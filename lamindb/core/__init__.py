"""Core library.

Settings & context:

.. autosummary::
   :toctree: .

   Context

Data loaders:

.. autosummary::
   :toctree: .

   MappedCollection

Modules:

.. autosummary::
   :toctree: .

   datasets
   storage
   logger

"""

from lamin_utils import logger
from lamin_utils._inspect import InspectResult

from .. import errors as exceptions
from . import datasets, types
from ._context import Context
from ._mapped_collection import MappedCollection

# backward compat
from ._settings import Settings
