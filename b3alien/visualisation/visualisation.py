import ee
import geemap
import folium
import tempfile
import webbrowser

def add_ee_layer(self, ee_object, vis_params, name):
    try:
        if isinstance(ee_object, ee.image.Image):
            map_id_dict = ee_object.getMapId(vis_params)
            folium.raster_layers.TileLayer(
                tiles=map_id_dict['tile_fetcher'].url_format,
                attr='Google Earth Engine',
                name=name,
                overlay=True,
                control=True
            ).add_to(self)
        elif isinstance(ee_object, ee.featurecollection.FeatureCollection):
            # Style the FeatureCollection (example with color)
            styled_fc = ee_object.style(**(vis_params or {'color': 'FF0000'}))
            map_id_dict = styled_fc.getMapId({})
            folium.raster_layers.TileLayer(
                tiles=map_id_dict['tile_fetcher'].url_format,
                attr='Google Earth Engine',
                name=name,
                overlay=True,
                control=True
            ).add_to(self)
        else:
            raise TypeError(f"Unsupported ee object type: {type(ee_object)}")
    except Exception as e:
        print("Could not add layer:", e)



# Patch Folium to support EE layers
def patch_folium():
    folium.Map.add_ee_layer = add_ee_layer


def visualize_ee_layers(layers, center=[0, 0], zoom=2, save_path=None, show=True):
    """
    Visualize multiple Earth Engine layers on a folium map.

    Parameters:
        layers: List of (ee_object, vis_params, layer_name)
        center: Center [lat, lon] for map
        zoom: Zoom level
        save_path: Optional filepath to save the HTML map
        show: If True, opens the map in the browser
    """
    m = folium.Map(location=center, zoom_start=zoom)

    for ee_object, vis_params, name in layers:
        m.add_ee_layer(ee_object, vis_params, name)

    folium.LayerControl().add_to(m)

    if save_path:
        m.save(save_path)
        print(f"Map saved to {save_path}")
        if show:
            webbrowser.open(f"file://{os.path.abspath(save_path)}")
    else:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as f:
            m.save(f.name)
            if show:
                webbrowser.open(f.name)