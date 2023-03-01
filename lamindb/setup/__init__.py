import lndb as _lndb
from lndb import *
from lndb import settings

from . import dev

__doc__ = _lndb.__doc__.replace("lndb", "lamindb.setup")
settings.__doc__ = settings.__doc__.replace("lamin", "lamindb.setup")
