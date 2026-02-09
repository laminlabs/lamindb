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
from ..examples import datasets  # backward compat
from . import subsettings, types
from ._context import Context
from ._settings import Settings


def __getattr__(name: str):
    # need to lazy import a few auxliary modules to maintain
    # backward compat
    # none of them should have been eagerly imported in the first place
    if name == "loaders":
        import importlib

        loaders = importlib.import_module(".loaders", package=__name__)
        return loaders
    if name == "storage":
        import importlib

        storage = importlib.import_module(".storage", package=__name__)
        return storage
    if name == "MappedCollection":
        from ._mapped_collection import MappedCollection

        return MappedCollection
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
