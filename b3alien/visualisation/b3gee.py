import ee
import json

def initialize(project):
    try:
        ee.Initialize()
    except Exception as e:
        ee.Authenticate()
        ee.Initialize(project=project)

def gdf_to_ee_featurecollection(gdf):
    geojson = json.loads(gdf.to_json())
    return ee.FeatureCollection(geojson)


