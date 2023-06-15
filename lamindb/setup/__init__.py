import lamindb_setup as _lamindb_setup
from lamindb_setup import *
from lamindb_setup import settings

from . import dev

__doc__ = _lamindb_setup.__doc__.replace("lamindb_setup", "lamindb.setup")
settings.__doc__ = settings.__doc__.replace("lamindb_setup", "lamindb.setup")
