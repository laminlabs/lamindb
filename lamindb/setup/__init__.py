import sys as _sys

import lamindb_setup as _lamindb_setup
from lamindb_setup import *  # noqa: F403
from lamindb_setup import (
    connect,
    delete,
    evolve,
    init,
    merge,
    settings,
    switch,
)

from . import core, errors, types

# Support the dotted import as well as access through ln.setup.evolve.
_sys.modules[f"{__name__}.evolve"] = evolve

del connect  # we have this at the root level, hence, we don't want it here
__doc__ = _lamindb_setup.__doc__.replace("lamindb_setup", "lamindb.setup")
settings.__doc__ = settings.__doc__.replace("lamindb_setup", "lamindb.setup")
