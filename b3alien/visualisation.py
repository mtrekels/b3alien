import folium
import tempfile
import webbrowser

# Patch Folium to support EE layers
def patch_folium():
    def add_ee_layer(self, ee_object, vis_params, name):
        map_id_dict = ee.Image(ee_object).getMapId(vis_params)
        folium.raster_layers.TileLayer(
            tiles=map_id_dict['tile_fetcher'].url_format,
            attr='Google Earth Engine',
            name=name,
            overlay=True,
            control=True
        ).add_to(self)
    folium.Map.add_ee_layer = add_ee_layer

def visualize_ee_layer(ee_object, vis_params, center=[0, 0], zoom=2, layer_name="Layer"):
    patch_folium()
    m = folium.Map(location=center, zoom_start=zoom)
    m.add_ee_layer(ee_object, vis_params, layer_name)
    folium.LayerControl().add_to(m)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as f:
        m.save(f.name)
        webbrowser.open(f.name)
