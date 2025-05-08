# tests/test_cube.py
import pytest
import numpy as np
import pandas as pd
from b3alien import b3cube as b3

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

        print(val)

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

