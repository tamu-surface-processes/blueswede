from . import gridding
from . import visualization
from . import preprocessing

# make available at top level like ANUGA does
from .preprocessing import read_polygon

# version stuff
from ._version import _version

__version__: str = _version()
