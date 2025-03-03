"""Core library.

CatManager:

.. autosummary::
   :toctree: .

   CatManager
   DataFrameCatManager
   AnnDataCatManager
   MuDataCatManager
   TiledbsomaCatManager
   CurateLookup

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

   loaders
   datasets
   storage
   logger

"""

from lamin_utils import logger
from lamin_utils._inspect import InspectResult

# from lamindb.models.query_manager import QueryManager
# from lamindb.models.query_set import QuerySet, RecordList
from lamindb.curators import (
    AnnDataCatManager,
    CatManager,
    CurateLookup,
    Curator,
    DataFrameCatManager,
    MuDataCatManager,
    TiledbsomaCatManager,
)

# below is backward compat
from lamindb.models import (
    BasicRecord,
    CanCurate,
    FeatureValue,
    HasParents,
    IsVersioned,
    ParamValue,
    Record,
    Registry,
    TracksRun,
    TracksUpdates,
    ValidateFields,
)

from .. import errors as exceptions
from . import datasets, fields, loaders, types
from ._context import Context
from ._mapped_collection import MappedCollection

# backward compat
from ._settings import Settings
