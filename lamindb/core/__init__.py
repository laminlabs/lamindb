"""Core library.

Settings & context:

.. autosummary::
   :toctree: .

   Settings
   subsettings
   Context

Artifact loaders:

.. autosummary::
   :toctree: .

   loaders

Data loaders:

.. autosummary::
   :toctree: .

   MappedCollection

Modules:

.. autosummary::
   :toctree: .

   storage
   logger

"""

from lamin_utils import logger
from lamin_utils._inspect import InspectResult

from .. import errors as exceptions  # backward compat
from ..base import types  # backward compat
from ..examples import datasets  # backward compat
from . import subsettings
from ._context import Context
from ._settings import Settings


def __getattr__(name: str):
    # need to lazy import a few auxliary modules to maintain backward compatibility
    # none of them should have been eagerly imported in the first place
    import importlib

    if name == "loaders":
        loaders = importlib.import_module(".loaders", package=__name__)
        globals()[name] = loaders
        return loaders
    if name == "storage":
        storage = importlib.import_module(".storage", package=__name__)
        globals()[name] = storage
        return storage
    if name == "MappedCollection":
        from ._mapped_collection import MappedCollection

        globals()[name] = MappedCollection
        return MappedCollection
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
