import lamindb_setup as _lamindb_setup
from lamindb_setup import *  # noqa: F403
from lamindb_setup import (
    connect,
    delete,
    init,
    settings,
)

from . import core

del connect  # we have this at the root level, hence, we don't want it here
__doc__ = _lamindb_setup.__doc__.replace("lamindb_setup", "lamindb.setup")
settings.__doc__ = settings.__doc__.replace("lamindb_setup", "lamindb.setup")
