import geopandas as gpd

def load_parquet(path):
    gdf = gpd.read_parquet(path)
    if gdf.empty:
        raise ValueError("GeoParquet file is empty or invalid.")
    return gdf
