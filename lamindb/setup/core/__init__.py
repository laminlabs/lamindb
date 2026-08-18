import lamindb_setup as _lamindb_setup
from lamindb_setup.core import *  # noqa: F403

__doc__ = _lamindb_setup.core.__doc__.replace("lamindb_setup", "lamindb.setup")


def __getattr__(name: str):
    return getattr(_lamindb_setup.core, name)
