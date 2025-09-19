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
import numpy as np


# Import the function and class from your module
from b3alien.b3cube import plot_richness, OccurrenceCube, find_correlations

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


class OccurrenceCubeStub:
    """Minimal stub that mimics the parts of OccurrenceCube your function needs."""
    def __init__(self, dataarray: xr.DataArray | None, df: pd.DataFrame):
        self.data = dataarray
        self.df = df


@pytest.fixture
def monkeypatch_xr_compute(monkeypatch):
    """
    Ensure xr.DataArray has a .compute() that returns an object whose .data has .todense().
    This mirrors the function's expectation when using dask/sparse arrays.
    """
    def _install():
        original_compute = getattr(xr.DataArray, "compute", None)

        def _dense_proxy(arr):
            class DenseProxy:
                def __init__(self, a): self._a = a
                def todense(self): return self._a
            return DenseProxy(arr)

        def fake_compute(self):
            # Return an object with dims/coords and .data.todense()
            class Dummy:
                def __init__(self, da):
                    self.dims = da.dims
                    self.coords = da.coords
                    # Use values to ensure a NumPy ndarray
                    self.data = _dense_proxy(da.values)
            return Dummy(self)

        # If .compute doesn't exist or is different, temporarily patch it
        monkeypatch.setattr(xr.DataArray, "compute", fake_compute, raising=False)
        return original_compute

    return _install


def test_find_correlations_basic_order_and_names(monkeypatch_xr_compute):
    # Install the safe compute() shim
    monkeypatch_xr_compute()

    # species keys and (partial) name map
    species_keys = [101, 102, 103]
    species_names = pd.DataFrame({
        "specieskey": [101, 102],           # intentionally omit 103 to test "Unknown (...)"
        "species": ["Oak", "Pine"],
    })

    # Build a tiny cube: dims (time=2, cell=2, species=3)
    # Presence (>0) layout:
    # t0,c0: 101,102
    # t0,c1: 101,103
    # t1,c0: 101,102
    # t1,c1: 103
    arr = np.zeros((2, 2, 3), dtype=int)
    # t0,c0 -> 101,102
    arr[0, 0, species_keys.index(101)] = 1
    arr[0, 0, species_keys.index(102)] = 1
    # t0,c1 -> 101,103
    arr[0, 1, species_keys.index(101)] = 1
    arr[0, 1, species_keys.index(103)] = 1
    # t1, c0 -> 101,102
    arr[1, 0, species_keys.index(101)] = 1
    arr[1, 0, species_keys.index(102)] = 1
    # t1, c1 -> 103
    arr[1, 1, species_keys.index(103)] = 1

    da = xr.DataArray(
        arr,
        dims=("time", "cell", "species"),
        coords={
            "time": [0, 1],
            "cell": [0, 1],
            "species": species_keys,
        },
        name="occurrence",
    )

    cube = OccurrenceCubeStub(da, species_names)

    result = find_correlations(cube, top_n=2)

    # Expected co-occurrence counts:
    # (101,102) -> 2  ; (101,103) -> 1  ; (102,103) -> 0
    # With top_n=2, we should get [(Oak, Pine), 2] and [(Oak, Unknown(103)), 1]
    assert len(result) == 2
    assert result[0] == (("Oak", "Pine"), 2)
    # Species 103 has no name in map -> "Unknown (103)"
    assert result[1] == (("Oak", "Unknown (103)"), 1)


def test_find_correlations_no_data_returns_empty(monkeypatch_xr_compute):
    monkeypatch_xr_compute()
    cube = OccurrenceCubeStub(dataarray=None, df=pd.DataFrame(columns=["specieskey", "species"]))
    assert find_correlations(cube) == []


def test_find_correlations_no_pairs_above_zero(monkeypatch_xr_compute):
    monkeypatch_xr_compute()
    # Each (time,cell) has only one species present -> no pair co-occurrences
    species_keys = [1, 2, 3, 4]
    arr = np.zeros((2, 2, 4), dtype=int)
    arr[0, 0, 0] = 1  # t0,c0 -> sp1
    arr[0, 1, 1] = 1  # t0,c1 -> sp2
    arr[1, 0, 2] = 1  # t1,c0 -> sp3
    arr[1, 1, 3] = 1  # t1,c1 -> sp4

    da = xr.DataArray(
        arr,
        dims=("time", "cell", "species"),
        coords={"time": [0, 1], "cell": [0, 1], "species": species_keys},
    )
    names = pd.DataFrame({"specieskey": species_keys, "species": [f"S{i}" for i in species_keys]})
    cube = OccurrenceCubeStub(da, names)

    assert find_correlations(cube) == []


def test_find_correlations_respects_top_n(monkeypatch_xr_compute):
    monkeypatch_xr_compute()
    # Construct cube with a clear ranking:
    # Pairs: (A,B)=3, (A,C)=2, (B,C)=1
    species_keys = [10, 20, 30]  # A, B, C
    names = pd.DataFrame({"specieskey": species_keys, "species": ["A", "B", "C"]})

    # Manually create (time,cell) blocks to realize those counts:
    # Use 3 blocks where A&B co-occur, 2 blocks where A&C co-occur, 1 block where B&C co-occur
    blocks = [
        {10, 20}, {10, 20}, {10, 20},  # three A&B
        {10, 30}, {10, 30},            # two A&C
        {20, 30},                      # one B&C
    ]
    T = len(blocks)
    C = 1
    S = len(species_keys)
    arr = np.zeros((T, C, S), dtype=int)
    for t, pair in enumerate(blocks):
        for key in pair:
            arr[t, 0, species_keys.index(key)] = 1

    da = xr.DataArray(
        arr,
        dims=("time", "cell", "species"),
        coords={"time": list(range(T)), "cell": [0], "species": species_keys},
    )
    cube = OccurrenceCubeStub(da, names)

    out_top2 = find_correlations(cube, top_n=2)
    assert out_top2 == [(("A", "B"), 3), (("A", "C"), 2)]

    out_top1 = find_correlations(cube, top_n=1)
    assert out_top1 == [(("A", "B"), 3)]