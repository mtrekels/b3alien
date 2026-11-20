import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely

#: CRS of the EEA reference grid (ETRS89-extended / LAEA Europe)
EEA_CRS = "EPSG:3035"

#: Regular expression for EEA reference grid cell codes, e.g. '1kmE2639N1762'
EEA_CELLCODE_PATTERN = r"^(?P<size>\d+)(?P<unit>m|km)E(?P<east>\d+)N(?P<north>\d+)$"

#: Column names used by GBIF cubes for the different grid systems
GRID_COLUMNS = (
    "cellcode",
    "eeacellcode",
    "eqdgccellcode",
    "eqdcellcode",
    "mgrscellcode",
    "isea3hcellcode",
)


def parse_eea_cellcode(codes):
    """
        Parse EEA reference grid cell codes into their resolution and lower-left corner.

        The EEA reference grid code has the form ``<size><unit>E<east>N<north>``,
        e.g. ``1kmE2639N1762``. The easting and northing are the coordinates of
        the lower-left corner in EPSG:3035, divided by 10^n, with n the number of
        trailing zeros of the cell size in metres (e.g. 1000 for 1km, 100 for 100m).

        Parameters
        ----------
        codes : array-like of str
            EEA cell codes.

        Returns
        -------
        pandas.DataFrame
            DataFrame (same index as ``codes`` if it is a Series) with the columns
            'cellcode', 'resolution' (cell size in metres), 'x' and 'y'
            (lower-left corner in metres, EPSG:3035).

        Raises
        ------
        ValueError
            If one or more codes are not valid EEA reference grid codes.
    """
    codes = pd.Series(codes, copy=False).astype(str)
    parts = codes.str.extract(EEA_CELLCODE_PATTERN)

    invalid = parts["size"].isna()
    if invalid.any():
        examples = codes[invalid].unique()[:5].tolist()
        raise ValueError(
            f"{int(invalid.sum())} value(s) are not valid EEA grid cell codes, e.g. {examples}"
        )

    size = parts["size"].astype("int64")
    resolution = np.where(parts["unit"] == "km", size * 1000, size).astype("int64")

    # Scaling factor 10^n, n = number of trailing zeros of the resolution in metres
    res_str = pd.Series(resolution).astype(str)
    trailing_zeros = (res_str.str.len() - res_str.str.rstrip("0").str.len()).to_numpy()
    factor = np.power(10, trailing_zeros).astype("int64")

    x = parts["east"].astype("int64").to_numpy() * factor
    y = parts["north"].astype("int64").to_numpy() * factor

    misaligned = (x % resolution != 0) | (y % resolution != 0)
    if misaligned.any():
        warnings.warn(
            f"{int(misaligned.sum())} EEA cell code(s) do not align with their grid resolution.",
            stacklevel=2,
        )

    return pd.DataFrame(
        {"cellcode": codes.to_numpy(), "resolution": resolution, "x": x, "y": y},
        index=codes.index,
    )


def is_eea_cellcode(codes):
    """
        Check whether all (non-null) values are EEA reference grid cell codes.

        Parameters
        ----------
        codes : array-like of str
            Cell codes to test.

        Returns
        -------
        bool
            True when every non-null value is a valid EEA cell code.
    """
    codes = pd.Series(codes, copy=False).dropna()
    if codes.empty:
        return False
    return bool(codes.astype(str).str.fullmatch(EEA_CELLCODE_PATTERN).all())


def eea_cellcode_to_geometry(codes):
    """
        Convert EEA reference grid cell codes into square polygons.

        Parameters
        ----------
        codes : array-like of str
            EEA cell codes (e.g. '10kmE263N176', '1kmE2639N1762', '100mE26390N17620').

        Returns
        -------
        geopandas.GeoSeries
            The grid cells as polygons in EPSG:3035.
    """
    parsed = parse_eea_cellcode(codes)
    res = parsed["resolution"].to_numpy()
    x = parsed["x"].to_numpy()
    y = parsed["y"].to_numpy()
    polygons = shapely.box(x, y, x + res, y + res)
    return gpd.GeoSeries(polygons, index=parsed.index, crs=EEA_CRS)


def eea_grid(codes):
    """
        Build a GeoDataFrame of the (unique) EEA grid cells.

        Parameters
        ----------
        codes : array-like of str
            EEA cell codes, duplicates are removed.

        Returns
        -------
        geopandas.GeoDataFrame
            One row per cell with the columns 'cellcode', 'resolution' and
            'geometry' (EPSG:3035).
    """
    unique_codes = pd.Series(pd.unique(pd.Series(codes).dropna()), dtype=str)
    parsed = parse_eea_cellcode(unique_codes)
    geometry = eea_cellcode_to_geometry(unique_codes)
    return gpd.GeoDataFrame(
        parsed[["cellcode", "resolution"]],
        geometry=geometry.values,
        crs=EEA_CRS,
    )


def add_eea_geometry(df, cell_col="eeacellcode"):
    """
        Add the EEA grid cell polygons to a (GBIF cube) dataframe.

        The geometry is computed once per unique cell and then mapped to all rows,
        which keeps this fast for large cubes.

        Parameters
        ----------
        df : pandas.DataFrame
            Dataframe containing a column with EEA cell codes.
        cell_col : str, optional
            Name of the column with EEA cell codes. Default is 'eeacellcode'.

        Returns
        -------
        geopandas.GeoDataFrame
            Copy of ``df`` with a 'geometry' column in EPSG:3035.
    """
    if cell_col not in df.columns:
        raise KeyError(f"Column '{cell_col}' not found in dataframe.")

    cells = pd.Categorical(df[cell_col].astype(str))
    cell_geoms = eea_cellcode_to_geometry(pd.Series(cells.categories)).to_numpy()
    geometry = cell_geoms[cells.codes]

    return gpd.GeoDataFrame(df, geometry=geometry, crs=EEA_CRS)


def to_geoparquet(csvFile, geoFile=None, crs="EPSG:4258", leftID='eqdcellcode', rightID='cellcode',
                  exportPath='./data/export.parquet', separator='\t'):
    """
        Convert a GBIF cube download into a GeoParquet file.

        The geometry is either taken from a GeoPackage (``geoFile``), or, when no
        GeoPackage is given and the cube uses the EEA reference grid, computed
        directly from the EEA cell codes (in EPSG:3035).

        Parameters
        ----------
        csvFile : str
            Path to the GBIF cube csv file (can be zipped).
        geoFile : str, optional
            Path to the GeoPackage file. If None, the geometry is derived from
            the EEA cell codes in column ``leftID``.
        crs : str, optional
            Coordinate reference system to be used for the GeoParquet file if the
            GeoPackage has none. Default is 'EPSG:4258'. Ignored for computed EEA
            geometries, which are always EPSG:3035.
        leftID : str, optional
            Column name within the GBIF cube to match the geometry. Default is 'eqdcellcode'.
            Use 'eeacellcode' for EEA grid cubes.
        rightID : str, optional
            Column name within the GeoPackage geometry. Default is 'cellcode'
        exportPath : str, optional
            Path to which the GeoParquet file needs to be exported.
        separator : str, optional
            Separator of the csv file. Default is tab.

        Returns
        -------
        A GeoParquet file at the location of exportPath
    """

    data = pd.read_csv(csvFile, sep=separator)

    if geoFile is None:
        if leftID not in data.columns:
            raise KeyError(
                f"Column '{leftID}' not found in the cube. Set leftID to the column with EEA cell codes."
            )
        gdf = add_eea_geometry(data, cell_col=leftID)
    else:
        geoRef = gpd.read_file(geoFile, engine='pyogrio', use_arrow=True)

        test_merge = pd.merge(data, geoRef, left_on=leftID, right_on=rightID)

        gdf = gpd.GeoDataFrame(test_merge, geometry='geometry')
        if gdf.crs is None:
            gdf.set_crs(crs, inplace=True)

    gdf.to_parquet(exportPath, engine="pyarrow", index=False)
