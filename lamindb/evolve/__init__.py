"""Public access to metadata schema evolution, implemented in lamindb-setup."""

import lamindb_setup.evolve as _evolve
from lamindb_setup.evolve import *  # noqa: F403
from lamindb_setup.evolve import __all__ as __all__
from lamindb_setup.evolve import __version__ as __version__

__doc__ = _evolve.__doc__
