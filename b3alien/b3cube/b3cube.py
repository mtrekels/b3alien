import os
import geopandas as gpd
import pandas as pd
import xarray as xr
import numpy as np
import sparse
import dask.array as da
import matplotlib.pyplot as plt
import matplotlib
import gcsfs
import folium
from folium import Choropleth
from shapely import wkt
from IPython.display import display
from b3alien.utils.runtime import in_jupyter
from b3alien.utils.geo import GRID_COLUMNS, EEA_CRS, is_eea_cellcode, add_eea_geometry
matplotlib.use("Agg")

#: Columns of a GBIF cube that are read as categories to reduce memory use
_CATEGORICAL_COLUMNS = ("kingdom", "phylum", "class", "order", "family", "genus",
                        "species", "yearmonth", "year") + GRID_COLUMNS

"""

Create a class object of the occurrence cube

"""


class OccurrenceCube():

    """
        Load a GeoParquet file (local or from GCS) into a sparse xarray cube.

        Parameters
        ----------
        filepath : str
            Path to the GeoParquet file (e.g. 'gs://bucket/file.parquet').
        dims : list or tuple, optional
            Dimension names. Default is ['time', 'cell', 'species'].
        coords : dict, optional
            Optional coordinates to assign to the cube.
        index_col : str or list, optional
            Column(s) to use for reshaping if needed.
        separator : str, optional
            Separator for GBIF CSV files. Default is tab ('\t').
        grid_column : str, optional
            Column holding the grid cell codes. By default this is detected
            automatically ('cellcode', 'eeacellcode', 'eqdgccellcode', ...).
            The column is renamed to 'cellcode' in self.df.
        add_geometry : bool, optional
            When loading a GBIF CSV cube on the EEA reference grid, compute the
            cell polygons (EPSG:3035) from the cell codes. Default is True.

        Returns
        -------
        b3cube.OccurrenceCube
            A sparse data cube loaded from the GeoParquet file.
            self.df contains a geopandas.DataFrame
            self.data a sparse xarray.Xarray
            self.grid_type is 'eea' for EEA reference grid cubes, otherwise None
            self.crs is the CRS of the cell geometries (None without geometry)

        Notes
        -----
        GBIF cubes on the new (Catalogue of Life based) backbone use
        alphanumeric keys (e.g. 'CPZZL'). Key columns are kept as strings in
        that case; cubes with purely numeric keys keep integer keys.
    """

    def __init__(self, filepath: str, source='geoparquet', gproject='', dims=None, coords=None, index_col=None,
                 separator='\t', grid_column=None, add_geometry=True):

        self.filepath = filepath
        self.source = source
        self.gproject = gproject
        self.dims = dims or ("time", "cell", "species")
        self.coords = coords
        self.index_col = index_col
        self.separator = separator
        self.grid_column = grid_column
        self.grid_type = None

        # Load Data
        if source == 'gbif':
            df = self._load_gbifcsv(filepath, separator)
        else:
            df = self._load_geoparquet(filepath, gproject)

        self.df = self._standardise_grid(df, grid_column, add_geometry)
        self.crs = getattr(self.df, "crs", None) if "geometry" in self.df.columns else None

        # Create cube
        self.data = self._create_xcube(self.df)

    def _load_gbifcsv(self, path, separator):
        """
        Load a GBIF CSV file (optionally zipped) from local disk using Pandas.

        Taxonomic names, time and grid columns are read as categories. Key
        columns are read as strings and only converted to integers when all
        keys are numeric (old GBIF backbone).
        """
        header = pd.read_csv(path, sep=separator, nrows=0).columns
        dtypes = {}
        for col in header:
            low = col.lower()
            if low.endswith("key") or low in _CATEGORICAL_COLUMNS:
                dtypes[col] = "category"

        df = pd.read_csv(path, sep=separator, dtype=dtypes)
        df.columns = df.columns.str.lower()

        for col in df.columns:
            if col.endswith("key") and isinstance(df[col].dtype, pd.CategoricalDtype):
                categories = df[col].cat.categories.astype(str)
                if len(categories) and categories.str.fullmatch(r"-?\d+").all():
                    df[col] = df[col].cat.rename_categories(categories.astype("int64"))

        return df

    def _standardise_grid(self, df, grid_column=None, add_geometry=True):
        """
        Make sure the cell codes are available in a 'cellcode' column and, for
        EEA reference grid cubes without geometry, add the cell polygons.
        """
        df.columns = df.columns.str.lower()

        if grid_column is not None:
            grid_column = grid_column.lower()
            if grid_column not in df.columns:
                raise KeyError(f"Grid column '{grid_column}' not found in the cube.")
        elif "cellcode" in df.columns:
            grid_column = "cellcode"
        else:
            grid_column = next((c for c in GRID_COLUMNS if c in df.columns), None)
            if grid_column is None:
                raise KeyError(
                    f"No grid cell column found. Expected one of {list(GRID_COLUMNS)}, "
                    "or pass grid_column explicitly."
                )

        if grid_column != "cellcode":
            if "cellcode" in df.columns:
                df = df.drop(columns="cellcode")
            df = df.rename(columns={grid_column: "cellcode"})
        self.grid_column = grid_column

        cells = df["cellcode"]
        unique_cells = cells.cat.categories if isinstance(cells.dtype, pd.CategoricalDtype) else cells.unique()
        if is_eea_cellcode(unique_cells):
            self.grid_type = "eea"
            if add_geometry and "geometry" not in df.columns:
                df = add_eea_geometry(df, cell_col="cellcode")

        return df

    def _load_geoparquet(self, path, gproject):
        """
        Load a GeoParquet file from local disk or GCS using GeoPandas. 
        """
        if path.startswith("gs://"):
            fs = gcsfs.GCSFileSystem(project=gproject)
            with fs.open(path) as f:
                gdf = gpd.read_parquet(f)
        else:
            gdf = gpd.read_parquet(path)

        if 'geometry' not in gdf.columns:
            raise ValueError("The input file must contain a 'geometry' column.")

        # Ensure geometry is parsed and valid
        gdf["geometry"] = gdf["geometry"].apply(wkt.loads) if gdf["geometry"].dtype == object else gdf["geometry"]
        gdf = gdf.set_geometry("geometry")

        return gdf

    def _create_xcube(self, df, dims=("time", "cell", "species")):
        """
        This function converts a GeoDataFrame into a sparse xarray cube with geometry metadata in case a GeoParquet file was loaded.
        In case of a pure GBIF cube, the geometry is ignored.  
        """
        # Convert to categorical
        df["yearmonth"] = pd.Categorical(df["yearmonth"])
        df["cellcode"] = pd.Categorical(df["cellcode"])
        df["specieskey"] = pd.Categorical(df["specieskey"])

        # Encode to integers
        time_codes = df["yearmonth"].cat.codes.values
        cell_codes = df["cellcode"].cat.codes.values
        species_codes = df["specieskey"].cat.codes.values

        # Build sparse cube
        sparse_cube = sparse.COO(
            coords=[time_codes, cell_codes, species_codes],
            data=df["occurrences"].astype("float32").values,
            shape=(
                df["yearmonth"].cat.categories.size,
                df["cellcode"].cat.categories.size,
                df["specieskey"].cat.categories.size
            )
        )

        # Prepare coordinates for the DataArray
        coords = {
            dims[0]: df["yearmonth"].cat.categories,
            dims[1]: df["cellcode"].cat.categories,
            dims[2]: df["specieskey"].cat.categories,
        }

        # Conditionally add geometry if the column exists
        if "geometry" in df.columns:
            cell_categories = df["cellcode"].cat.categories
            geometry_per_cell = df.drop_duplicates("cellcode").set_index("cellcode").loc[cell_categories]["geometry"]
            coords["geometry"] = (dims[1], geometry_per_cell.values)

        # Create xarray DataArray
        cube = xr.DataArray(
            sparse_cube,
            dims=self.dims,
            coords=coords,
            name="occurrences"
        )

        return cube

    def _species_richness(self, normalized=False):
        """
            Calculate species richness per cell from the sparse cube.

            Parameters
            ----------
            normalized : bool
                Whether to compute normalized richness (richness / total occurrences).

            Returns
            -------
            None
                Sets self.richness as a DataFrame with columns 'cell' and 'richness' or 'normalized_richness'.
        """
        # 1. Binary presence
        presence = (self.data > 0)

        # 2. Collapse time dimension using logical OR → was the species *ever* seen in this cell?
        presence_any_time = presence.any(dim="time")  # shape: (cell, species)

        # 3. Sum species per cell (species richness)
        species_richness = presence_any_time.sum(dim="species")  # shape: (cell,)

        total_occurrences = self.data.sum(dim=["time", "species"])

        if not normalized:
            # 4. Get the non-zero values and indices
            coords = species_richness.data.coords  # (1D arrays of indices)
            values = species_richness.data.data    # the richness values

            # 5. Convert integer cell indices to real labels (from .coords['cell'])
            cell_labels = species_richness.coords["cell"].values

            richness_df = pd.DataFrame({
                "cell": cell_labels[coords[0]],
                "richness": values
            })

            self.richness = richness_df

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

            self.richness = norm_df

    def _filter_species(self, speciesKey):
        """
            Filter the cube to only include data for a specific speciesKey.

            Parameters
            ----------
            speciesKey : int or str
                The GBIF speciesKey to filter on.

            Returns
            -------
            None
                Updates self.df and self.data to only include the specified speciesKey.
        """

        self.df = self.df[self.df['specieskey'].eq(speciesKey)]
        self.data = self.data.sel(species=speciesKey)


def find_correlations(cube: OccurrenceCube, top_n=10):
    """
    Calculates species co-occurrence using efficient matrix operations.

    This method identifies pairs of species that occur in the same 'cell'
    at the same 'time'.

    Parameters
    ----------
        cube (OccurrenceCube): An instance of the OccurrenceCube class.
        top_n (int): The number of top correlated pairs to return.

    Returns
    -------
        list: A list of tuples, where each tuple contains the species names
              and their co-occurrence count, sorted in descending order.
    """
    if cube.data is None:
        print("Data cube is not loaded. Cannot find correlations.")
        return []
        
    # 1. Create a boolean cube for presence (occurrence > 0)
    # Using int32 to prevent overflow when summing large arrays.
    presence = (cube.data > 0).astype('int32')

    # 2. Compute the co-occurrence matrix using broadcasting.
    # This multiplies the presence array with itself along a new 'species_copy' dimension
    # and sums over the 'time' and 'cell' dimensions.
    # The result is a (species x species) matrix.
    co_occurrence = (presence * presence.rename({'species': 'species_copy'})).sum(dim=['time', 'cell'])

    # 3. If the data is a Dask array, trigger the computation.
    # Then, explicitly densify the sparse array before converting to pandas.
    computed_co = co_occurrence.compute()
    df_co = xr.DataArray(
        computed_co.data.todense(),
        dims=computed_co.dims,
        coords=computed_co.coords
    ).to_pandas()
    df_co.columns.name = 'species'
    df_co.index.name = 'species'

    # 4. Extract the upper triangle to avoid duplicate pairs (A,B vs B,A)
    # and self-pairs (A,A). k=1 excludes the diagonal.
    df_upper = df_co.where(np.triu(np.ones(df_co.shape), k=1).astype(bool))
    
    # 5. Stack the matrix to get a Series of (species, species) -> count
    stacked_pairs = df_upper.stack().sort_values(ascending=False)

    # Filter out pairs with zero co-occurrences
    strongest_pairs = stacked_pairs[stacked_pairs > 0]
    
    if strongest_pairs.empty:
        print("No co-occurring species pairs found in the provided dataset.")
        return []

    # 6. Create a mapping from specieskey to species name for readability
    species_map = cube.df.drop_duplicates('specieskey').set_index('specieskey')['species']
    
    # 7. Format the output
    results = []
    for (key1, key2), count in strongest_pairs.head(top_n).items():
        name1 = species_map.get(key1, f"Unknown ({key1})")
        name2 = species_map.get(key2, f"Unknown ({key2})")
        results.append(((name1, name2), int(count)))
        
    return results


def filter_time_window(df, start_year, end_year, cols=["year", "rate"]):
    """Filter a two column dataframe based on year window.
    
    Parameters
    ----------
    df: pandas.DataFrame
        A two column dataframe with one time axs
    start_year : np.int
        Integer number of the earliest date to filter on.
    end_year : np.int
        Integer number of the latest date to filter on.
    cols: list of str, optional
        List of the names of the columns of the dataframe
    

    Returns
    -------
    time: pandas.Series
        Filtered time series
    filtered_var: pandas.Series
        Filtered variable corresponding to the time series
    
    """
    filtered = df[(df[cols[0]] >= start_year) & (df[cols[0]] <= end_year)].reset_index(drop=True)
    return filtered[cols[0]], filtered[cols[1]]


def aggregate_count_per_cell(cube, taxonRank, taxon):
    """
    Aggregate the counts per taxonomic level per cell. This can be used as a
    normalization factor

    Parameters
    ----------
    cube: OccurrenceCube
        An OccurrenceCube with geometries added
    taxonRank: dwc.taxonRank term
        Level at which the aggregation needs to be performed
    taxon: str
        Name of the taxon at which aggregation needs to be performed
    plot: bool, optional
        whether the aggregated count per cell needs to be plotted on a map

    Returns
    -------
    gdf: GeoDataFrame
        dataframe containing geometries and the aggregated occurrence count

    """
    # Step 1: Filter by the desired class
    filtered_df = cube.df[cube.df[taxonRank] == taxon]

    # Step 2: Select only the relevant columns
    count_str = taxonRank + "count"
    columns_to_keep = ["yearmonth", "cellcode", taxonRank, count_str, "geometry"]
    subset_df = filtered_df[columns_to_keep]

    # Step 3: Drop duplicate rows based on those columns
    gdf = subset_df.drop_duplicates()

    return gdf


def plot_richness(cube, normalized=False, html_path='richness_map.html', default_crs="EPSG:4326"):
    """
    Plot species richness from an OccurrenceCube instance, handling CRS dynamically.

    Parameters
    ----------
    cube : OccurrenceCube
        An instance of the OccurrenceCube class with a loaded data cube and optional geometry.
    normalized : bool, optional
        Whether to plot normalized richness (richness / total occurrences). Default is False.
    html_path : str, optional
        File path to save the HTML map if not in Jupyter. Default is 'richness
        _map.html'.
    default_crs : str, optional
        The default CRS to assign if none is detected. Default is "EPSG:4326".

    Returns
    -------
    folium.Map or geemap.Map
        A map object visualizing species richness, either for display in Jupyter or saved as HTML for external viewing. 
    """
    # 1. Ensure richness is computed
    if not hasattr(cube, "richness"):
        cube._species_richness(normalized=normalized)

    # 2. Extract richness and geometry
    richness_df = cube.richness.copy()
    geometry = cube.data.coords["geometry"].values
    cell_labels = cube.data.coords["cell"].values
    
    # Create the mapping for geometry alignment
    idx_map = dict(zip(cell_labels, range(len(cell_labels))))
    
    # 3. Create GeoDataFrame
    gdf_plot = gpd.GeoDataFrame(
        richness_df,
        geometry=geometry[richness_df["cell"].map(idx_map)],
        crs=_cube_crs(cube)  # None if unknown; default_crs is assigned below
    )

    # 4. Handle CRS logic
    # If the cube already has a CRS defined, use it; otherwise, assign the default.
    if gdf_plot.crs is None:
        print(f"No CRS detected. Setting to default: {default_crs}")
        gdf_plot.set_crs(default_crs, inplace=True)
    
    # Folium REQUIRES EPSG:4326. If it's something else, we reproject.
    if gdf_plot.crs != "EPSG:4326":
        print(f"Reprojecting from {gdf_plot.crs} to EPSG:4326 for Folium compatibility.")
        gdf_plot = gdf_plot.to_crs("EPSG:4326")

    # 5. Optionally simplify geometry (helps with performance)
    gdf_plot["geometry"] = gdf_plot["geometry"].simplify(0.01)

    # 6. Center map using the reprojected coordinates
    # centroid = gdf_plot.geometry.union_all().centroid # GeoPandas 1.0+
    # For older versions, use: gdf_plot.geometry.unary_union.centroid
    total_bounds = gdf_plot.total_bounds # [minx, miny, maxx, maxy]
    center_lat = (total_bounds[1] + total_bounds[3]) / 2
    center_lon = (total_bounds[0] + total_bounds[2]) / 2
    
    m = folium.Map(location=[center_lat, center_lon], zoom_start=7)

    # 7. Plotting
    richness_col = 'richness' if not normalized else 'normalized_richness'
    
    Choropleth(
        geo_data=gdf_plot,
        name='Species Richness',
        data=gdf_plot,
        columns=['cell', richness_col],
        key_on='feature.properties.cell',
        fill_color='YlGnBu',
        fill_opacity=0.7,
        line_opacity=0.2,
        legend_name='Species Richness'
    ).add_to(m)

    folium.LayerControl().add_to(m)

    # Output handling
    if 'in_jupyter' in globals() and in_jupyter():
        return m
    else:
        m.save(html_path)
        print(f"Map saved as HTML to {os.path.abspath(html_path)}")


def _align_keys(cube_keys, species_to_keep):
    """
        Cast checklist keys to the type of the species keys in the cube, so that
        numeric (old backbone) and alphanumeric (new backbone) keys compare correctly.
    """
    cube_keys = np.asarray(cube_keys)
    if np.issubdtype(cube_keys.dtype, np.integer):
        aligned = []
        for key in species_to_keep:
            try:
                aligned.append(int(key))
            except (TypeError, ValueError):
                continue
        return aligned
    return [str(key) for key in species_to_keep]


def _cube_crs(cube):
    """Return the CRS of the cell geometries of a cube, or None if unknown."""
    crs = getattr(cube, "crs", None)
    if crs is None:
        crs = getattr(getattr(cube, "df", None), "crs", None)
    if crs is None and getattr(cube, "grid_type", None) == "eea":
        crs = EEA_CRS
    return crs


def grid_geodataframe(cube, column=None, aggfunc="sum"):
    """
    Build a GeoDataFrame with one row per grid cell of an OccurrenceCube.

    Parameters
    ----------
    cube : OccurrenceCube
        Cube with cell geometries (GeoParquet input, or a GBIF CSV cube on the
        EEA reference grid).
    column : str, optional
        Value to attach to each cell:

        - 'richness' or 'normalized_richness': species richness per cell;
        - 'occurrences_total': total number of occurrences per cell;
        - any other column of ``cube.df`` (e.g. 'distinctobservers',
          'kingdomcount'), aggregated per cell with ``aggfunc``.

        If None, only the cell geometries are returned.
    aggfunc : str or callable, optional
        Aggregation applied to a ``cube.df`` column per cell. Default is 'sum'.

    Returns
    -------
    geopandas.GeoDataFrame
        Columns 'cell', the optional value column and 'geometry'.
    """
    if "geometry" not in cube.data.coords:
        raise ValueError(
            "The cube has no cell geometries. Load a GeoParquet cube, or a GBIF cube on the EEA grid."
        )

    gdf = gpd.GeoDataFrame(
        {"cell": np.asarray(cube.data.coords["cell"].values)},
        geometry=list(cube.data.coords["geometry"].values),
        crs=_cube_crs(cube),
    )

    if column is None:
        return gdf

    if column in ("richness", "normalized_richness"):
        cube._species_richness(normalized=(column == "normalized_richness"))
        values = cube.richness.set_index("cell")[column]
    elif column == "occurrences_total":
        totals = cube.data.sum(dim=["time", "species"])
        values = pd.Series(np.asarray(totals.data.todense() if hasattr(totals.data, "todense") else totals.data),
                           index=np.asarray(totals.coords["cell"].values))
    elif column in cube.df.columns:
        values = cube.df.groupby("cellcode", observed=True)[column].agg(aggfunc)
    else:
        raise KeyError(f"Column '{column}' not found in the cube.")

    gdf[column] = gdf["cell"].map(values)
    return gdf


def _to_grid_gdf(source, column=None, aggfunc="sum"):
    """Convert the input of plot_grid to a GeoDataFrame."""
    if isinstance(source, gpd.GeoDataFrame):
        return source
    if isinstance(source, OccurrenceCube) or hasattr(source, "data"):
        return grid_geodataframe(source, column=column, aggfunc=aggfunc)

    # Assume an iterable of EEA cell codes
    from b3alien.utils.geo import eea_grid
    gdf = eea_grid(source).rename(columns={"cellcode": "cell"})
    return gdf


def plot_grid(source, column=None, kind="static", log=False, cmap="YlGnBu", aggfunc="sum",
              title=None, legend_label=None, ax=None, output_path=None, html_path="grid_map.html",
              basemap=False, edgecolor=None, linewidth=0.05, figsize=(10, 12)):
    """
    Visualise a grid (e.g. the EEA reference grid), optionally coloured by a value per cell.

    Works on GBIF cubes on the EEA grid without any GeoParquet or GeoPackage:
    the cell polygons are derived from the EEA cell codes.

    Parameters
    ----------
    source : OccurrenceCube, geopandas.GeoDataFrame or iterable of str
        The cube to plot, a GeoDataFrame with cell geometries, or a list of EEA
        cell codes (e.g. ``cube.df['cellcode'].unique()``).
    column : str, optional
        Value to colour the cells with. For a cube this can be 'richness',
        'normalized_richness', 'occurrences_total' or any column of cube.df
        (aggregated per cell with ``aggfunc``). For a GeoDataFrame, a column
        of that dataframe. If None, only the grid cells are drawn.
    kind : {'static', 'interactive'}, optional
        'static' draws a matplotlib map in the native CRS (EPSG:3035 for the
        EEA grid); 'interactive' returns a folium map. Default is 'static'.
    log : bool, optional
        Use a logarithmic colour scale (static maps only). Default is False.
    cmap : str, optional
        Colour map (matplotlib for static, ColorBrewer name for interactive).
    aggfunc : str or callable, optional
        Aggregation per cell for cube.df columns. Default is 'sum'.
    title : str, optional
        Title of the map.
    legend_label : str, optional
        Label of the colour bar / legend. Defaults to the column name.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on (static maps only).
    output_path : str, optional
        If given, the static figure is saved to this path.
    html_path : str, optional
        Path to save the interactive map when not running in Jupyter.
    basemap : bool, optional
        Add a basemap with contextily, if installed (static maps only).
    edgecolor : str, optional
        Colour of the cell borders. Default is no borders when colouring by
        ``column`` (best for fine grids) and 'grey' when drawing the grid only.
    linewidth : float, optional
        Width of the cell borders. Default is 0.05.
    figsize : tuple, optional
        Figure size of the static map.

    Returns
    -------
    (matplotlib.figure.Figure, matplotlib.axes.Axes) or folium.Map
        The static figure and axes, or the interactive folium map.
    """
    gdf = _to_grid_gdf(source, column=column, aggfunc=aggfunc)
    if gdf.crs is None:
        raise ValueError("The grid has no CRS defined. Set it with gdf.set_crs(...).")
    if column is not None and column not in gdf.columns:
        raise KeyError(f"Column '{column}' not found.")

    label = legend_label or column

    if kind == "interactive":
        return _plot_grid_interactive(gdf, column, cmap, label, html_path)
    if kind != "static":
        raise ValueError("kind must be 'static' or 'interactive'.")

    import matplotlib.colors as mcolors

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    if column is None:
        gdf.plot(ax=ax, facecolor="none", edgecolor=edgecolor or "grey", linewidth=max(linewidth, 0.2), zorder=2)
    else:
        plot_gdf = gdf.copy()
        plot_gdf[column] = pd.to_numeric(plot_gdf[column], errors="coerce")
        norm = None
        if log:
            positive = plot_gdf.loc[plot_gdf[column] > 0, column]
            if positive.empty:
                raise ValueError(f"Column '{column}' has no positive values for a log scale.")
            plot_gdf = plot_gdf[plot_gdf[column] > 0]
            norm = mcolors.LogNorm(vmin=positive.min(), vmax=positive.max())
        plot_gdf.plot(
            column=column, cmap=cmap, norm=norm, ax=ax, legend=True,
            edgecolor=edgecolor or "none", linewidth=linewidth if edgecolor else 0, zorder=2,
            legend_kwds={"label": label, "shrink": 0.6},
            missing_kwds={"color": "lightgrey"},
        )

    if basemap:
        try:
            import contextily as cx
            cx.add_basemap(ax, crs=gdf.crs.to_string(), source=cx.providers.CartoDB.Positron, zorder=1)
        except Exception as e:  # contextily missing or tiles unavailable
            print(f"Basemap could not be added: {e}")

    ax.set_aspect("equal")
    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=14)
    ax.text(0.99, 0.01, f"CRS: {gdf.crs.to_string()}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color="grey")

    if output_path:
        fig.savefig(output_path, bbox_inches="tight", dpi=300)

    return fig, ax


def _plot_grid_interactive(gdf, column, cmap, label, html_path):
    """Interactive folium version of plot_grid."""
    if len(gdf) > 50000:
        print(f"Warning: {len(gdf)} cells; the interactive map may be slow. "
              "Consider kind='static'.")

    gdf_4326 = gdf.to_crs("EPSG:4326")
    bounds = gdf_4326.total_bounds
    m = folium.Map(location=[(bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2], zoom_start=7)
    m.fit_bounds([[bounds[1], bounds[0]], [bounds[3], bounds[2]]])

    if column is None:
        folium.GeoJson(
            gdf_4326[["cell", "geometry"]],
            name="Grid",
            style_function=lambda _: {"color": "#3366cc", "weight": 0.5, "fillOpacity": 0.05},
            tooltip=folium.GeoJsonTooltip(fields=["cell"]),
        ).add_to(m)
    else:
        data = gdf_4326[["cell", column, "geometry"]].dropna(subset=[column])
        Choropleth(
            geo_data=data,
            name=label,
            data=data,
            columns=["cell", column],
            key_on="feature.properties.cell",
            fill_color=cmap,
            fill_opacity=0.7,
            line_opacity=0.2,
            legend_name=label,
        ).add_to(m)

    folium.LayerControl().add_to(m)

    if not in_jupyter():
        m.save(html_path)
        print(f"Map saved as HTML to {os.path.abspath(html_path)}")
    return m


def cumulative_species(cube, species_to_keep):

    """
        Calculate the cumulative number of species in a OccurrenceCube.

        Parameters
        ----------
        cube : b3alien.b3cube.OccurrenceCube
            Species OccurrenceCube from GBIF.
        species_to_keep : numpy.array
            Array of GBIF speciesKeys that need to be taken into account to calculate the cumulative species number of a cube.
        geom : str, optional

        Returns
        -------
        df1 : pandas.DataFrame
            Sparse dataframe that still contains the cumulative sum per grid cell.
        df2 : pandas.DataFrame
            Cumulative dataframe cell independent.
    """

    # Wrap sparse array in Dask array with one or more chunks
    dask_sparse_array = da.from_array(cube.data.data, chunks=(100, 100, 1000))  # tune chunking for your use case

    # Replace data in cube
    cube_dask_sparse = cube.data.copy(data=dask_sparse_array)

    species_mask = cube_dask_sparse["species"].isin(
        _align_keys(cube_dask_sparse["species"].values, species_to_keep)
    )
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

    df_time["time"] = pd.to_datetime(df_time["time"], format="%Y-%m", errors="coerce")

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
    df_cumulative["time"] = pd.to_datetime(df_cumulative["time"], format="%Y-%m", errors="coerce")

    return df_sparse, df_cumulative


def plot_cumsum(df_cumulative):
    """
        Create a plot of the cumulative number of species.

        Parameters
        ----------
        df_cumulative : pandas.DataFrame
            Datagrame containing the cumulative number over tume.

        Returns
        -------
        matplotlib.plot
            A plot of the cumulative number of species.
    """

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
    """
        Only count a species established when it is present in more than one cell.

        Parameters
        ----------
        df_sparse : pandas.DataFrame
            Datagrame containing the species richness per grid cell.

        Returns
        -------
        pandas.DataFrame
            Cumulative species when in multiple cells.
    """
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

    df_cumulative_cells = pd.DataFrame(cumulative,
                                       columns=["time",
                                                "cumulative_species_2cells"])

    return df_cumulative_cells

def filter_multiple_occ(df_sparse):
    """
        Only count a species established when there are multiple occurrences in a cell.

        Parameters
        ----------
        df_sparse : pandas.DataFrame
            Datagrame containing the species richness per grid cell.

        Returns
        -------
        pandas.DataFrame
            Cumulative species when multiple occurrences in a cell.
    """
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

    df_cumulative_occ = pd.DataFrame(cumulative, 
                                     columns=["time",
                                              "cumulative_species_2occ"])

    return df_cumulative_occ

def calculate_rate(df_cumulative):
    """
        Calculate the rate of establishment from the cumulative distribution.

        Parameters
        ----------
        df_cumulative : pandas.DataFrame
            Datagrame containing the cumulative distribution.

        Returns
        -------
        s1 : pandas.Series
            Series of the time axis.
        s2 : pandas.Series
            Series of the rate of establishment.
    """
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

        # Both cumulative values are taken at the end of their year, so the
        # elapsed time is the number of years between them (> 1 if years are missing).
        species_diff = current_year_data['cumulative_species'] - previous_year_data['cumulative_species']
        time_diff_years = current_year_data.name - previous_year_data.name

        annual_rate = species_diff / time_diff_years
        annual_rate_gbif.append(annual_rate)
        annual_time_gbif.append(current_year_data.name)  # Year is the index after groupby

    annual_time_gbif = [int(year) for year in annual_time_gbif]

    # Convert the lists to Pandas Series for easier plotting
    annual_rate_gbif_series = pd.Series(annual_rate_gbif)
    annual_time_gbif_series = pd.Series(annual_time_gbif)

    return annual_time_gbif, annual_rate_gbif


def get_survey_effort(cube, dateFormat='%Y-%m', calc_type='total'):
    """
        Estimate the survey effort in an OccurrenceCube.

        Parameters
        ----------
        cube : b3alien.b3cube.OccurrenceCube
            Species OccurrenceCube from GBIF.
        dateFormat : str, optional
            Dateformat stored in the OccurrenceCube. Default is '%Y-%m'
        calc_type : str, optional
            Type of survey effort to be calculated. 
                'distinct' : total number of distinct observers
                'total' : total number of occurrences
                Default is total.

        Returns
        -------
        df : pandas.DataFrame
            Dataframe containing time and the chosen measurement for survey effort.
    """
    if calc_type == 'distinct':
        # Group by 'yearmonth' and sum
        distinct_observers_over_time = cube.df.groupby('yearmonth', observed=True)['distinctobservers'].sum()

        # Convert the index to datetime (if it's not already)
        distinct_observers_over_time.index = pd.to_datetime(
            distinct_observers_over_time.index, format=dateFormat, errors='coerce'
        )

        # Resample to yearly frequency using new 'YE' standard
        distinct_observers_yearly = distinct_observers_over_time.resample('YE').sum()

        # Filter for years from 1900 onward
        distinct_observers_filtered = distinct_observers_yearly[distinct_observers_yearly.index.year >= 1900]

        # Convert to DataFrame
        df = distinct_observers_filtered.reset_index()
        df.columns = ['date', 'distinct_observers']  # Rename columns for clarity

        return df

    else:

        total_occurrences_over_time = cube.data.sum(dim=['cell', 'species'])

        # Convert the time coordinates to datetime objects if they aren't already
        total_occurrences_over_time['time'] = pd.to_datetime(total_occurrences_over_time['time'].values, format=dateFormat, errors='coerce')
        df = pd.DataFrame({
            'time': total_occurrences_over_time['time'].values,
            'total_occurrences': total_occurrences_over_time.data.data
        })

        # (Optional) Drop rows with invalid or missing time
        df = df.dropna(subset=['time'])

        return df

