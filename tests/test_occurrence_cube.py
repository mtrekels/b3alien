# test_occurrence_cube.py

import os
import pytest
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import Polygon
import folium
from unittest.mock import MagicMock
import sparse

# Import the function and class from your module
from b3alien.b3cube import plot_richness, OccurrenceCube

@pytest.fixture
def mock_occurrence_cube():
    """
    Creates a mock OccurrenceCube object with the necessary attributes 
    (.data, .richness) for testing the plot_richness function.
    """
    # 1. Define some simple geometries and cell labels
    poly1 = Polygon([(0, 0), (1, 1), (1, 0)])
    poly2 = Polygon([(2, 2), (3, 3), (3, 2)])
    cell_labels = ['cell_A', 'cell_B']
    geometries = [poly1, poly2]

    # 2. Create the .richness dataframe that the function expects
    richness_df = pd.DataFrame({
        "cell": ["cell_A", "cell_B"],
        "richness": [10, 25]
    })
    sparse_cube = sparse.COO(
        coords=[[0, 0], [0, 1], [0, 0]], # time, cell, species indices
        data=[10, 25],
        shape=(1, 2, 1) # time, cell, species shape
    )
    # 3. Create the .data xarray.DataArray with geometry in its coordinates
    data_array = xr.DataArray(
        sparse_cube,  # Dummy data values
        dims=("time", "cell", "species"),
        coords={
            "time": ['2020-01'],
            "cell": cell_labels,
            "species": [1001],
            "geometry": ("cell", geometries)
        }
    )

    # 4. Use MagicMock to create a fake OccurrenceCube instance
    # The 'spec=OccurrenceCube' ensures it mimics the real class
    mock_cube = MagicMock(spec=OccurrenceCube)
    mock_cube.data = data_array
    mock_cube.richness = richness_df
    mock_cube._species_richness = MagicMock() # Also mock the internal method

    return mock_cube

# test_occurrence_cube.py (continued)

def test_plot_richness_in_jupyter(mocker, mock_occurrence_cube):
    """
    Verifies that plot_richness calls display() when in a Jupyter environment.
    """
    # CORRECTED PATH: Patch the function where it is used.
    mock_in_jupyter = mocker.patch('b3alien.b3cube.b3cube.in_jupyter', return_value=True)
    
    mock_display = mocker.patch('b3alien.b3cube.b3cube.display')
    mock_save = mocker.patch('folium.Map.save')

    # --- Call the function ---
    plot_richness(mock_occurrence_cube)

    # --- Assertions ---
    mock_in_jupyter.assert_called_once()
    mock_display.assert_called_once()
    mock_save.assert_not_called()

def test_plot_richness_as_script(mocker, mock_occurrence_cube):
    """
    Verifies that plot_richness calls m.save() when not in a Jupyter environment.
    """
    # CORRECTED PATH: Patch the function where it is used.
    mock_in_jupyter = mocker.patch('b3alien.b3cube.b3cube.in_jupyter', return_value=False)
    
    mock_display = mocker.patch('b3alien.b3cube.b3cube.display')
    mock_save = mocker.patch('folium.Map.save')

    # --- Call the function ---
    test_path = 'my_test_map.html'
    plot_richness(mock_occurrence_cube, html_path=test_path)

    # --- Assertions ---
    mock_in_jupyter.assert_called_once()
    mock_save.assert_called_once_with(test_path)
    mock_display.assert_not_called()

## Test 3: (Bonus) Verifies it calculates richness if missing
def test_plot_richness_recreates_richness_correctly(mocker, mock_occurrence_cube):
    """
    Verifies that plot_richness, when .richness is missing, calls the
    real _species_richness method and correctly calculates the result.
    """
    # Setup:
    # 1. Attach the REAL method to our mock cube instance.
    #    This binds the method to the instance, so 'self' will work correctly.
    mock_occurrence_cube._species_richness = OccurrenceCube._species_richness.__get__(mock_occurrence_cube)

    # 2. Delete the attribute to trigger the calculation
    del mock_occurrence_cube.richness
    
    # Mock the plotting functions since we only care about the calculation
    mocker.patch('b3alien.b3cube.b3cube.in_jupyter', return_value=True)
    mocker.patch('b3alien.b3cube.b3cube.display')

    # --- Call the function that triggers the calculation ---
    plot_richness(mock_occurrence_cube)

    # --- Assertions ---
    # Now, instead of checking if a mock was called, we check the actual result!
    assert hasattr(mock_occurrence_cube, 'richness')
    
    # Check the calculated values. Our sparse data had one species in each
    # of the two cells, so richness for both should be 1.
    expected_df = pd.DataFrame({
        "cell": ["cell_A", "cell_B"],
        "richness": [1, 1]
    })
    pd.testing.assert_frame_equal(
        mock_occurrence_cube.richness.reset_index(drop=True),
        expected_df.reset_index(drop=True)
    )