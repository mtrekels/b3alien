import geopandas as gpd
import pandas as pd
import xarray as xr
import sparse
import dask.array as da
import numpy as np
import matplotlib as plt

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

def plot_richness(richness_df, gdf_from_gcs, geom='eqdgccellcode'):

    gdf_plot = pd.merge(richness_df, gdf_from_gcs, left_on='cell', right_on=geom)

    gdf_plot = gpd.GeoDataFrame(gdf_plot, geometry="geometry", crs=gdf_from_gcs.crs)

    fig, ax = plt.subplots(figsize=(10, 10))
    gdf_plot.plot(
        column="richness",
        cmap="viridis",
        legend=True,
        linewidth=0.1,
        edgecolor="grey",
        ax=ax
    )
    ax.set_title("Species Richness per QDGC Cell")
    ax.axis("off")
    plt.show()


def cumulative_species(cube, species_to_keep):
    # Wrap sparse array in Dask array with one or more chunks
    dask_sparse_array = da.from_array(cube.data, chunks=(100, 100, 1000))  # tune chunking for your use case

    # Replace data in cube
    cube_dask_sparse = cube.copy(data=dask_sparse_array)

    species_mask = cube_dask_sparse["species"].isin(species_to_keep)
    filtered_cube = cube_dask_sparse.where(species_mask, drop=True)

    # Grab the underlying sparse.COO object from Dask
    sparse_block = filtered_cube.data.compute()  # Warning: loads full filtered cube into RAM!

    # Extract sparse coordinates
    coords = sparse_block.coords  # shape: (ndim, nnz)
    data = sparse_block.data      # non-zero values

    # Map indices to labels
    time_labels = filtered_cube.coords["time"].values
    species_labels = filtered_cube.coords["species"].values
    cell_labels = filtered_cube.coords["cell"].values

    # Use the sparse indices to create a DataFrame
    df_sparse = pd.DataFrame({
        "time": time_labels[coords[0]],
        "cell": cell_labels[coords[1]],
        "species": species_labels[coords[2]],
        "occurrences": sparse_block.data
    })

    # Drop duplicates and compute cumulative species count
    df_sparse = df_sparse.drop_duplicates()
    df_sparse["seen"] = 1
    df_time = (
        df_sparse.groupby("time")["species"]
        .nunique()
        .cumsum()
        .reset_index(name="cumulative_species")
    )

    df_time["time"] = pd.to_datetime(df_time["time"], format="%Y-%M", errors="coerce")

    # fix to have the real cumsum
    # Step 1: Remove duplicates (species × time)
    df_sparse_unique = df_sparse[["time", "species"]].drop_duplicates()

    # Step 2: Sort by time
    df_sparse_unique = df_sparse_unique.sort_values("time")

    # Step 3: Track cumulative species using a set
    seen_species = set()
    cumulative = []

    for time, group in df_sparse_unique.groupby("time"):
        new_species = set(group["species"])
        seen_species.update(new_species)
        cumulative.append((time, len(seen_species)))

    # Step 4: Create cumulative DataFrame
    df_cumulative = pd.DataFrame(cumulative, columns=["time", "cumulative_species"])
    df_cumulative["time"] = pd.to_datetime(df_cumulative["time"], format="%Y-%M", errors="coerce")

    return df_sparse, df_cumulative

def plot_cumsum(df_cumulative):
    
    plt.figure(figsize=(10, 5))
    plt.plot(df_cumulative["time"], df_cumulative["cumulative_species"], marker="o")
    plt.title("Cumulative Unique Species Over Time")
    plt.xlabel("Time")
    plt.ylabel("Unique Species Observed")
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def filter_multiple_cells(df_sparse):
    # Ensure cell is in your DataFrame
    assert "cell" in df_sparse.columns

    # 1. Count unique cells per (time, species)
    species_cell_counts = (
        df_sparse.groupby(["time", "species"])["cell"]
        .nunique()
        .reset_index(name="cell_count")
    )

    # 2. Keep only species seen in at least 2 cells
    species_multi_cell = species_cell_counts[species_cell_counts["cell_count"] >= 2]

    # 3. Track cumulative set of species across time
    species_multi_cell = species_multi_cell.sort_values("time")

    seen_species = set()
    cumulative = []

    for time, group in species_multi_cell.groupby("time"):
        new_species = set(group["species"])
        seen_species.update(new_species)
        cumulative.append((time, len(seen_species)))

    df_cumulative_cells = pd.DataFrame(cumulative, columns=["time", "cumulative_species_2cells"])

    return df_cumulative_cells

def filter_multiple_occ(df_sparse):
    # Ensure 'occurrences' and 'cell' are present
    assert "occurrences" in df_sparse.columns and "cell" in df_sparse.columns

    # 1. Total occurrences per (time, species, cell)
    species_cell_occ = (
        df_sparse.groupby(["time", "species", "cell"])["occurrences"]
        .sum()
        .reset_index()
    )

    # 2. Filter for species that had ≥ 2 occurrences in any cell
    species_with_2occ = (
        species_cell_occ[species_cell_occ["occurrences"] >= 2]
        .drop_duplicates(subset=["time", "species"])
    )

    # 3. Cumulative species count logic
    species_with_2occ = species_with_2occ.sort_values("time")

    seen_species = set()
    cumulative = []

    for time, group in species_with_2occ.groupby("time"):
        new_species = set(group["species"])
        seen_species.update(new_species)
        cumulative.append((time, len(seen_species)))

    df_cumulative_occ = pd.DataFrame(cumulative, columns=["time", "cumulative_species_2occ"])

    return df_cumulative_occ

def calculate_rate(df_cumulative):
    # --- Processing GBIF data (Monthly) to get an approximate annual rate ---
    df_cumulative["time"] = pd.to_datetime(df_cumulative["time"])
    df_cumulative_rate = df_cumulative.sort_values(by="time").copy()

    # Group data by year and calculate the total species difference for each year
    annual_data = df_cumulative_rate.groupby(df_cumulative_rate['time'].dt.year).agg(
        cumulative_species=('cumulative_species', 'last'),  # Get the last cumulative species value for the year
        first_time=('time', 'first')  # Get the first timestamp for the year
    )

    # Calculate annual rate using the grouped data
    annual_rate_gbif = []
    annual_time_gbif = []

    for i in range(1, len(annual_data)):
        current_year_data = annual_data.iloc[i]
        previous_year_data = annual_data.iloc[i - 1]

        species_diff = current_year_data['cumulative_species'] - previous_year_data['cumulative_species']
        time_diff_years = (current_year_data['first_time'] - previous_year_data['first_time']).days / 365.25

        annual_rate = species_diff / time_diff_years
        annual_rate_gbif.append(annual_rate)
        annual_time_gbif.append(current_year_data.name)  # Year is the index after groupby

    annual_time_gbif = [int(year) for year in annual_time_gbif]

    # Convert the lists to Pandas Series for easier plotting
    annual_rate_gbif_series = pd.Series(annual_rate_gbif)
    annual_time_gbif_series = pd.Series(annual_time_gbif)

    return annual_time_gbif, annual_rate_gbif

'''
This function is not finished yet

def get_survey_effort(cube, calc_type='do'):

    if calc_type == 'do':
        # Assuming your DataFrame is named 'gdf_from_gcs'
        distinct_observers_over_time = gdf_from_gcs.groupby('yearmonth')['distinctobservers'].sum()

        # Convert the 'yearmonth' index to datetime objects if it isn't already
        distinct_observers_over_time.index = pd.to_datetime(distinct_observers_over_time.index, format='%Y-%m', errors='coerce')

        # Resample to yearly frequency and sum
        distinct_observers_yearly = distinct_observers_over_time.resample('Y').sum()

        # Filter data from 1940 onwards
        distinct_observers_filtered = distinct_observers_yearly[distinct_observers_yearly.index.year >= 1900]

        return distinct_observers_filtered

    else:

        total_occurrences_over_time = cube.sum(dim=['cell', 'species'])

        # Convert the time coordinates to datetime objects if they aren't already
        total_occurrences_over_time['time'] = pd.to_datetime(total_occurrences_over_time['time'].values, format='%Y-%m', errors='coerce')

        return total_occurrences_over_time
    
'''