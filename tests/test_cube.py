# tests/test_cube.py
import pytest
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from b3alien import b3cube as b3
from b3alien import griis
from b3alien import utils
import folium

def test_cube_loading():
    cube = b3.OccurrenceCube("tests/data/data_PT-30.parquet")
    assert cube.data.dims == ("time", "cell", "species")
    assert "geometry" in cube.data.coords


def test_cube_content():
    cube = b3.OccurrenceCube("tests/data/data_PT-30.parquet")
    data = cube.data
    # Check shape
    assert data.shape == (
        len(data.coords["time"]),
        len(data.coords["cell"]),
        len(data.coords["species"])
    )
    # Example: check a specific value
    # Replace with actual expected values from your test data
    expected_time = "2018-07"
    expected_cell = "W017N32BBDD"
    expected_species = 2979000
    expected_occurrences = 3.0

    if (expected_time in data.coords["time"].values and
        expected_cell in data.coords["cell"].values and
        expected_species in data.coords["species"].values):

        val = data.drop_vars("geometry").sel(
            time=expected_time,
            cell=expected_cell,
            species=expected_species
        ).item()

        assert np.isclose(val, expected_occurrences), f"Expected {expected_occurrences}, got {val}"

def test_richness():
    cube = b3.OccurrenceCube("tests/data/data_PT-30.parquet")
    
    cube._species_richness()

    assert isinstance(cube.richness, pd.DataFrame)

def test_survey_effort():

    cube = b3.OccurrenceCube("tests/data/data_PT-30.parquet")
    total = b3.get_survey_effort(cube)
    distinctObs = b3.get_survey_effort(cube, calc_type='distinct')

    assert (total['total_occurrences'] >= 0).all()
    assert (distinctObs['distinct_observers'] >= 0).all()

def test_counts_per_cell():
    cube = cube = b3.OccurrenceCube("tests/data/data_PT-30.parquet")
    gdf = b3.aggregate_count_per_cell(cube, "genus", "Oxalis")

    assert isinstance(gdf, gpd.GeoDataFrame)

def test_occurrence_cube_init():
    cube = b3.OccurrenceCube("tests/data/data_PT-30.parquet")
    assert hasattr(cube, "data")
    assert isinstance(cube.data, xr.DataArray)

def test_occurrence_cube_gdf():
    cube = b3.OccurrenceCube("tests/data/data_PT-30.parquet")
    assert hasattr(cube, "df")
    assert isinstance(cube.df, gpd.GeoDataFrame)

def test_cumulative_species():
    cube = b3.OccurrenceCube("tests/data/data_PT-30.parquet")
    checklist = griis.CheckList("tests/data/dwca-griis-portugal-v1.9/merged_distr.txt")
    df_sparse, df_cumulative = b3.cumulative_species(cube, checklist.species)
    assert isinstance(df_sparse, pd.DataFrame)
    assert isinstance(df_cumulative, pd.DataFrame)
    assert df_cumulative.shape[0] > 0

def test_calculate_rate():
    cube = b3.OccurrenceCube("tests/data/data_PT-30.parquet")
    checklist = griis.CheckList("tests/data/dwca-griis-portugal-v1.9/merged_distr.txt")
    df_sparse, df_cumulative = b3.cumulative_species(cube, checklist.species)
    annual_time, annual_rate = b3.calculate_rate(df_cumulative)
    assert isinstance(annual_time, list)
    assert isinstance(annual_rate, list)
    assert len(annual_time) == len(annual_rate)