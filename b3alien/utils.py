"""
Additional utilities

"""
import pandas as pd
import geopandas as gpd
import pyarrow


def detect_runtime():
    """
    Detects the runtime environment where the code is executed.

    Returns:
        str: One of "Jupyter Notebook", "IPython Terminal", or "Standard Python Script"
    """
    try:
        from IPython import get_ipython
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return "Jupyter Notebook"
        elif shell == 'TerminalInteractiveShell':
            return "IPython Terminal"
        else:
            return "Other IPython"
    except (ImportError, AttributeError, NameError):
        return "Standard Python Script"

def in_jupyter():
    """Returns True if running inside a Jupyter Notebook or Lab."""
    return detect_runtime() == "Jupyter Notebook"

def in_ipython():
    """Returns True if running inside any IPython shell (not standard Python)."""
    return detect_runtime() != "Standard Python Script"

def in_script():
    """Returns True if running in a standard Python script (non-interactive)."""
    return detect_runtime() == "Standard Python Script"

def to_geoparquet(csvFile, geoFile, leftID=eqdcellcode, rightID=cellCode, exportPath='./data/export.parquet'):
    """
       Convert a GBIF cube download into a GeoParquet file, using the geometry of a GPKG

        Parameters
        ----------
        csvFile : str
            Path to the GBIF cube csv file.
        geoFile : str
            Path to the GeoPackage file.
        leftID : str, optional
            Column name within the GBIF cube to match the geometry. Default is 'edqcellcode'.
        rightID : str, optional
            Column name within the GeoPackage geometry. Default is 'cellCode'
        exportPath : str, optional
            Path to which the GeoParquet file needs to be exported.

        Returns
        -------
        A GeoParquet file at the location of exportPath
    """

    data = pd.read_csv(csvFile, sep='\t')
    geoRef = gpd.read_file(geoFile, engine='pyogrio', use_arrow=True, crs="EPSG:4326")

    test_merge = pd.merge(data, qdgc_ref, left_on=leftID, right_on=rightID)

    gdf = gpd.GeoDataFrame(test_merge, geometry='geometry')
    if gdf.crs is None:
        gdf.set_crs(crs, inplace=True) 

    gdf.to_parquet(exportPath, engine="pyarrow", index=False)

