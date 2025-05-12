import pandas as pd
import geopandas as gpd
import pyarrow

def to_geoparquet(csvFile, geoFile, leftID='eqdcellcode', rightID='cellCode', exportPath='./data/export.parquet'):
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