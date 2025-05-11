import ee
import json

def initialize(project):
    """
        Inititialise a Google project

        Parameters
        ----------
        project : str
            Name of your Google project
    """
    try:
        ee.Initialize()
    except Exception as e:
        ee.Authenticate()
        ee.Initialize(project=project)

def gdf_to_ee_featurecollection(gdf):
    """
        Transform a GeoDataFrame in to an Earth Engine Feature collection

        Parameters
        ----------
        gdf : geopandas.DataFrame
            Species OccurrenceCube from GBIF.
            
        Returns
        -------
        ee.FeatureCollection
    """
    geojson = json.loads(gdf.to_json())
    return ee.FeatureCollection(geojson)


