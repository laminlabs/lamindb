from typing import TYPE_CHECKING

import lamindb_setup as _lamindb_setup
from lamindb_setup import *  # noqa: F403
from lamindb_setup import (
    connect,
    delete,
    init,
    merge,
    settings,
    switch,
)

from . import core, errors, types

if TYPE_CHECKING:
    from . import evolve as evolve

del connect  # we have this at the root level, hence, we don't want it here
__doc__ = _lamindb_setup.__doc__.replace("lamindb_setup", "lamindb.setup")
settings.__doc__ = settings.__doc__.replace("lamindb_setup", "lamindb.setup")


def __getattr__(name: str):
    if name == "evolve":
        from importlib import import_module

        return import_module(".evolve", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
