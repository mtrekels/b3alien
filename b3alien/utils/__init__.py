"""

Utils
=====

"""

from .runtime import detect_runtime
from .runtime import in_jupyter
from .runtime import in_ipython
from .runtime import in_script
from .geo import to_geoparquet

__all__ = ["detect_runtime", "in_jupyter", "in_ipython", "in_script", "to_geoparquet"]