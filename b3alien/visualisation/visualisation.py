import ee
import geemap
import folium
import tempfile
import webbrowser
import os

from b3alien.utils.runtime import in_jupyter


def add_ee_layer(self, ee_object, vis_params, name):
    """
        Add Earth Engine layers to a Folium map.
    """
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


def patch_folium():
    """
     Monkey-patch folium.Map to support EE layers.
    """
    if not hasattr(folium.Map, 'add_ee_layer'):
        folium.Map.add_ee_layer = add_ee_layer


def visualize_ee_layers(layers, center=[0, 0], zoom=2, save_path=None, show=True, image=False):
    """
    Visualize Earth Engine layers either in Jupyter (geemap) or save as HTML/PNG via folium.

    Parameters:
        layers (list): List of tuples (ee_object, vis_params, name)
        center (list): [lat, lon] center of map
        zoom (int): Zoom level
        save_path (str): Optional filepath to save HTML
        show (bool): Show in browser (script) or inline (Jupyter)
        image (bool): Attempt to save a PNG image (requires Selenium)

    Returns:
        folium.Map or geemap.Map
    """
    patch_folium()

    if in_jupyter():
        m = geemap.Map(center=center, zoom=zoom)
        for ee_object, vis_params, name in layers:
            m.addLayer(ee_object, vis_params, name)
        m.addLayerControl()
        if show:
            display(m)
        return m

    else:
        m = folium.Map(location=center, zoom_start=zoom)
        for ee_object, vis_params, name in layers:
            m.add_ee_layer(ee_object, vis_params, name)
        folium.LayerControl().add_to(m)

        if save_path is None:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as f:
                save_path = f.name

        m.save(save_path)
        print(f"Map saved to {save_path}")

        if show:
            webbrowser.open(f"file://{os.path.abspath(save_path)}")

        if image:
            try:
                from selenium import webdriver
                from selenium.webdriver.chrome.options import Options
                from PIL import Image
                import time

                options = Options()
                options.headless = True
                options.add_argument("--window-size=1200,800")
                driver = webdriver.Chrome(options=options)
                driver.get(f"file://{os.path.abspath(save_path)}")
                time.sleep(3)
                img_path = save_path.replace(".html", ".png")
                driver.save_screenshot(img_path)
                driver.quit()
                print(f"Screenshot saved to {img_path}")
            except Exception as e:
                print("Screenshot failed:", e)

        return m
