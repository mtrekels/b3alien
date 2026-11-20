# tests/test_eea_grid.py
import zipfile

import numpy as np
import pandas as pd
import geopandas as gpd
import pytest
import matplotlib
import folium

from b3alien import b3cube, griis, utils
from b3alien.b3cube import OccurrenceCube

matplotlib.use("Agg")


# ===========================================================================
# EEA cell code parsing / geometry
# ===========================================================================

@pytest.mark.parametrize("code, res, x, y", [
    ("1kmE2639N1762", 1000, 2639000, 1762000),
    ("10kmE263N176", 10000, 2630000, 1760000),
    ("100kmE26N17", 100000, 2600000, 1700000),
    ("100mE26390N17620", 100, 2639000, 1762000),
    ("25mE2639000N1762000", 25, 2639000, 1762000),
    ("250mE263900N176200", 250, 2639000, 1762000),
])
def test_parse_eea_cellcode(code, res, x, y):
    parsed = utils.parse_eea_cellcode([code]).iloc[0]
    assert parsed["resolution"] == res
    assert parsed["x"] == x
    assert parsed["y"] == y


def test_parse_eea_cellcode_invalid():
    with pytest.raises(ValueError):
        utils.parse_eea_cellcode(["1kmE2639N1762", "W017N32BBDD"])


def test_is_eea_cellcode():
    assert utils.is_eea_cellcode(["1kmE2639N1762", "10kmE263N176"])
    assert not utils.is_eea_cellcode(["W017N32BBDD"])
    assert not utils.is_eea_cellcode([])


def test_eea_cellcode_to_geometry():
    geoms = utils.eea_cellcode_to_geometry(["1kmE2639N1762", "10kmE263N176"])
    assert geoms.crs == utils.EEA_CRS
    assert geoms.iloc[0].bounds == (2639000, 1762000, 2640000, 1763000)
    assert np.isclose(geoms.iloc[1].area, 10000 ** 2)


def test_eea_grid_unique_cells():
    grid = utils.eea_grid(["1kmE1N1", "1kmE1N1", "1kmE1N2"])
    assert isinstance(grid, gpd.GeoDataFrame)
    assert len(grid) == 2
    assert set(grid.columns) == {"cellcode", "resolution", "geometry"}


def test_add_eea_geometry():
    df = pd.DataFrame({"eeacellcode": ["1kmE1N1", "1kmE1N2", "1kmE1N1"], "v": [1, 2, 3]})
    gdf = utils.add_eea_geometry(df)
    assert gdf.crs == utils.EEA_CRS
    assert gdf.geometry.iloc[0].equals(gdf.geometry.iloc[2])
    assert gdf.geometry.iloc[1].bounds == (1000, 2000, 2000, 3000)


# ===========================================================================
# OccurrenceCube on a GBIF CSV with EEA grid and new-backbone keys
# ===========================================================================

def write_eea_cube(tmp_path, zipped=False):
    df = pd.DataFrame({
        "kingdom": ["Animalia"] * 5,
        "kingdomkey": ["N"] * 5,
        "species": ["Ardea ibis", "Ardea ibis", "Ardea cinerea", "Ardea cinerea", "Ardea alba"],
        "specieskey": ["CPZZL", "CPZZL", "GCHS", "GCHS", "4321"],
        "yearmonth": ["2020-01", "2021-03", "2020-01", "2022-05", "2022-05"],
        "eeacellcode": ["1kmE2639N1762", "1kmE2640N1762", "1kmE2639N1762", "1kmE2641N1762", "1kmE2641N1762"],
        "kingdomcount": [5, 3, 5, 2, 2],
        "distinctobservers": [1, 2, 1, 3, 3],
        "occurrences": [2, 1, 4, 1, 1],
    })
    path = tmp_path / "cube.csv"
    df.to_csv(path, sep="\t", index=False)
    if zipped:
        zpath = tmp_path / "cube.zip"
        with zipfile.ZipFile(zpath, "w") as zf:
            zf.write(path, arcname="0000000-000000000000000.csv")
        return zpath
    return path


@pytest.mark.parametrize("zipped", [False, True])
def test_eea_cube_loading(tmp_path, zipped):
    cube = OccurrenceCube(str(write_eea_cube(tmp_path, zipped)), source="gbif")

    assert cube.grid_type == "eea"
    assert cube.grid_column == "eeacellcode"
    assert cube.crs == utils.EEA_CRS
    assert isinstance(cube.df, gpd.GeoDataFrame)
    assert "cellcode" in cube.df.columns
    assert "geometry" in cube.data.coords
    assert cube.data.shape == (3, 3, 3)  # 3 months x 3 cells x 3 species
    # Alphanumeric keys of the new backbone are kept as strings
    assert set(cube.data.coords["species"].values) == {"CPZZL", "GCHS", "4321"}


def test_eea_cube_without_geometry(tmp_path):
    cube = OccurrenceCube(str(write_eea_cube(tmp_path)), source="gbif", add_geometry=False)
    assert cube.grid_type == "eea"
    assert "geometry" not in cube.df.columns
    assert cube.crs is None


def test_numeric_keys_stay_numeric(tmp_path):
    df = pd.DataFrame({
        "specieskey": [2979000, 2979001],
        "yearmonth": ["2020-01", "2020-02"],
        "eeacellcode": ["1kmE1N1", "1kmE1N2"],
        "occurrences": [1, 2],
    })
    path = tmp_path / "cube.csv"
    df.to_csv(path, sep="\t", index=False)
    cube = OccurrenceCube(str(path), source="gbif")
    assert np.issubdtype(np.asarray(cube.data.coords["species"].values).dtype, np.integer)


def test_cumulative_species_with_string_keys(tmp_path):
    cube = OccurrenceCube(str(write_eea_cube(tmp_path)), source="gbif")
    # 4321 is given as int but stored as str in the cube
    d_s, d_c = b3cube.cumulative_species(cube, ["CPZZL", 4321])
    assert set(d_s["species"]) == {"CPZZL", "4321"}
    assert d_c["cumulative_species"].iloc[-1] == 2


def test_cumulative_species_parses_months(tmp_path):
    cube = OccurrenceCube(str(write_eea_cube(tmp_path)), source="gbif")
    _, d_c = b3cube.cumulative_species(cube, ["CPZZL", "GCHS", "4321"])
    expected = pd.to_datetime(["2020-01-01", "2021-03-01", "2022-05-01"])
    assert list(d_c["time"]) == list(expected)
    assert (d_c["time"].dt.minute == 0).all()


def test_eea_cube_richness(tmp_path):
    cube = OccurrenceCube(str(write_eea_cube(tmp_path)), source="gbif")
    gdf = b3cube.grid_geodataframe(cube, column="richness")
    richness = dict(zip(gdf["cell"], gdf["richness"]))
    assert richness == {"1kmE2639N1762": 2, "1kmE2640N1762": 1, "1kmE2641N1762": 2}


# ===========================================================================
# Grid visualisation
# ===========================================================================

def test_grid_geodataframe_columns(tmp_path):
    cube = OccurrenceCube(str(write_eea_cube(tmp_path)), source="gbif")

    grid = b3cube.grid_geodataframe(cube)
    assert len(grid) == 3 and grid.crs == utils.EEA_CRS

    obs = b3cube.grid_geodataframe(cube, column="distinctobservers")
    assert dict(zip(obs["cell"], obs["distinctobservers"]))["1kmE2641N1762"] == 6

    occ = b3cube.grid_geodataframe(cube, column="occurrences_total")
    assert occ["occurrences_total"].sum() == 9

    with pytest.raises(KeyError):
        b3cube.grid_geodataframe(cube, column="does_not_exist")


@pytest.mark.parametrize("column, log", [(None, False), ("richness", False), ("distinctobservers", True)])
def test_plot_grid_static(tmp_path, column, log):
    cube = OccurrenceCube(str(write_eea_cube(tmp_path)), source="gbif")
    out = tmp_path / "map.png"
    fig, ax = b3cube.plot_grid(cube, column=column, log=log, output_path=str(out))
    assert out.exists()
    matplotlib.pyplot.close(fig)


def test_plot_grid_from_cellcodes(tmp_path):
    fig, ax = b3cube.plot_grid(["10kmE263N176", "10kmE264N176"])
    assert len(ax.collections) == 1
    matplotlib.pyplot.close(fig)


def test_plot_grid_interactive(tmp_path, mocker):
    mocker.patch("b3alien.b3cube.b3cube.in_jupyter", return_value=False)
    cube = OccurrenceCube(str(write_eea_cube(tmp_path)), source="gbif")
    html = tmp_path / "map.html"
    m = b3cube.plot_grid(cube, column="richness", kind="interactive", html_path=str(html))
    assert isinstance(m, folium.Map)
    assert html.exists()


def test_plot_grid_invalid_kind(tmp_path):
    cube = OccurrenceCube(str(write_eea_cube(tmp_path)), source="gbif")
    with pytest.raises(ValueError):
        b3cube.plot_grid(cube, kind="3d")


def test_to_geoparquet_without_geopackage(tmp_path):
    out = tmp_path / "cube.parquet"
    utils.to_geoparquet(str(write_eea_cube(tmp_path)), leftID="eeacellcode", exportPath=str(out))
    gdf = gpd.read_parquet(out)
    assert gdf.crs == utils.EEA_CRS
    cube = OccurrenceCube(str(out))
    assert cube.grid_type == "eea"


# ===========================================================================
# ColDP checklists and the new backbone
# ===========================================================================

def write_coldp(tmp_path):
    folder = tmp_path / "coldp"
    folder.mkdir()
    pd.DataFrame({
        "col:ID": ["1", "2", "3"],
        "col:scientificName": ["Ardea ibis", "Ardea cinerea", "Unknownia fakea"],
        "col:authorship": ["Linnaeus, 1758", "", ""],
        "col:rank": ["species"] * 3,
        "col:kingdom": ["Animalia"] * 3,
    }).to_csv(folder / "NameUsage.tsv", sep="\t", index=False)
    pd.DataFrame({
        "col:taxonID": ["1", "2", "3"],
        "col:area": ["Portugal"] * 3,
        "col:establishmentMeans": ["introduced"] * 3,
    }).to_csv(folder / "Distribution.tsv", sep="\t", index=False)
    pd.DataFrame({
        "col:taxonID": ["1", "1", "2"],
        "col:property": ["dwc:habitat", "gbif:isInvasive", "dwc:habitat"],
        "col:value": ["terrestrial|freshwater", "true", "marine"],
    }).to_csv(folder / "TaxonProperty.tsv", sep="\t", index=False)
    return folder


def test_is_coldp(tmp_path):
    folder = write_coldp(tmp_path)
    assert griis.is_coldp(folder)
    zpath = tmp_path / "coldp.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        for f in folder.iterdir():
            zf.write(f, arcname=f.name)
    assert griis.is_coldp(zpath)
    assert not griis.is_coldp(tmp_path / "missing")


def test_do_taxon_matching_coldp(tmp_path, mocker):
    folder = write_coldp(tmp_path)
    lookup = {"Ardea ibis": ["CPZZL"], "Ardea cinerea": ["GCHS"]}
    mocker.patch(
        "b3alien.griis.griis._match_v2",
        side_effect=lambda name, *a, **k: lookup.get(name, ["Uncertain"]),
    )

    out = tmp_path / "out"
    merged = griis.do_taxon_matching(str(folder), output_dir=str(out))
    assert len(merged) == 3

    profile = pd.read_csv(out / "speciesprofile.txt", sep="\t")
    row1 = profile.set_index("id").loc[1]
    assert row1["isInvasive"] == "Invasive"
    assert row1["habitat"] == "Terrestrial|Freshwater"

    cl = griis.CheckList(str(out))
    assert sorted(cl.species) == ["CPZZL", "GCHS"]

    cl_inv = griis.CheckList(str(out), use_profile=True, invasive_val="Invasive")
    assert cl_inv.species == ["CPZZL"]


def test_match_v2_resolves_to_species(mocker):
    response = mocker.Mock()
    response.json.return_value = {
        "usage": {"key": "5JZ3T", "rank": "SUBSPECIES"},
        "classification": [{"key": "N", "rank": "KINGDOM"}, {"key": "493JT", "rank": "SPECIES"},
                           {"key": "5JZ3T", "rank": "SUBSPECIES"}],
        "diagnostics": {"matchType": "EXACT"},
    }
    mocker.patch("requests.get", return_value=response)
    assert griis.get_speciesKey("Olea europaea subsp. europaea",
                                checklistKey=griis.COL_XR_CHECKLIST_KEY) == ["493JT"]

    response.json.return_value = {"usage": {"key": "6B6T", "rank": "GENUS"},
                                  "classification": [], "diagnostics": {"matchType": "HIGHERRANK"}}
    assert griis.get_speciesKey("Oxalis foo", checklistKey=griis.COL_XR_CHECKLIST_KEY) == ["Uncertain"]
