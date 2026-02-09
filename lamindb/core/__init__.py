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

from .. import errors as exceptions
from ..examples import datasets  # backward compat
from . import loaders, subsettings, types
from ._context import Context
from ._settings import Settings


def __getattr__(name: str):
    # need to lazy import MappedCollection - it should likely moved out of the
    # root namespace over time
    if name == "MappedCollection":
        from ._mapped_collection import MappedCollection

        return MappedCollection
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
