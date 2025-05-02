import geopandas as gpd
import pandas as pd
import xarray as xr
import sparse
import dask.array as da
import numpy as np

def create_xcube(df):
    # Convert categorical dimensions
    df["yearmonth"] = pd.Categorical(df["yearmonth"])
    df["cellCode"] = pd.Categorical(df["cellCode"])
    df["specieskey"] = pd.Categorical(df["specieskey"])

    # Extract code arrays (integer labels) and levels
    time_codes = df["yearmonth"].cat.codes.values
    cell_codes = df["cellCode"].cat.codes.values
    species_codes = df["specieskey"].cat.codes.values

    # Sparse 3D cube (occurrence > 0)
    sparse_cube = sparse.COO(
        coords=[time_codes, cell_codes, species_codes],
        data=df["occurrences"].astype("float32").values,  # or just 1 for presence/absence
        shape=(
            df["yearmonth"].cat.categories.size,
            df["cellCode"].cat.categories.size,
            df["specieskey"].cat.categories.size
        )
    )

    # Wrap into Xarray
    cube = xr.DataArray(
        sparse_cube,
        dims=("time", "cell", "species"),
        coords={
            "time": df["yearmonth"].cat.categories,
            "cell": df["cellCode"].cat.categories,
            "species": df["specieskey"].cat.categories
        },
        name="occurrences"
    )

    return(cube)

def species_richness(cube, normalized=False):
    # 1. Binary presence
    presence = (cube > 0)

    # 2. Collapse time dimension using logical OR → was the species *ever* seen in this cell?
    presence_any_time = presence.any(dim="time")  # shape: (cell, species)

    # 3. Sum species per cell (species richness)
    species_richness = presence_any_time.sum(dim="species")  # shape: (cell,)

    total_occurrences = cube.sum(dim=["time", "species"])

    if normalized == False:
        # 4. Get the non-zero values and indices
        coords = species_richness.data.coords  # (1D arrays of indices)
        values = species_richness.data.data    # the richness values

        # 5. Convert integer cell indices to real labels (from .coords['cell'])
        cell_labels = species_richness.coords["cell"].values

        richness_df = pd.DataFrame({
            "cell": cell_labels[coords[0]],
            "richness": values
        })

        return richness_df

    else:
        epsilon = 1e-6
        normalized_richness = species_richness / (total_occurrences + epsilon)

        coords = normalized_richness.data.coords
        values = normalized_richness.data.data
        cell_labels = normalized_richness.coords["cell"].values

        # Build a DataFrame
        norm_df = pd.DataFrame({
            "cell": cell_labels[coords[0]],
            "normalized_richness": values
        })

        return norm_df



    