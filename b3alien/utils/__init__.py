"""

Utils
=====

"""

from .runtime import detect_runtime
from .runtime import in_jupyter
from .runtime import in_ipython
from .runtime import in_script
from .geo import to_geoparquet
from .geo import parse_eea_cellcode
from .geo import is_eea_cellcode
from .geo import eea_cellcode_to_geometry
from .geo import eea_grid
from .geo import add_eea_geometry
from .geo import EEA_CRS

__all__ = ["detect_runtime", "in_jupyter", "in_ipython", "in_script", "to_geoparquet",
           "parse_eea_cellcode", "is_eea_cellcode", "eea_cellcode_to_geometry",
           "eea_grid", "add_eea_geometry", "EEA_CRS"]
