"""Core library.

Registries:

.. autosummary::
   :toctree: .

   BasicRecord
   Record
   Registry
   QuerySet
   QueryManager
   RecordList
   FeatureManager
   ParamManager
   LabelManager
   IsVersioned
   CanCurate
   HasParents
   TracksRun
   TracksUpdates
   ParamValue
   FeatureValue
   InspectResult
   ValidateFields
   fields

Curators:

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

   Settings
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
   subsettings
   logger

"""

from lamin_utils import logger
from lamin_utils._inspect import InspectResult

# from lamindb.models.query_manager import QueryManager
# from lamindb.models.query_set import QuerySet, RecordList
# from lamindb.core._feature_manager import FeatureManager, ParamManager
# from lamindb.core._label_manager import LabelManager
# from lamindb.curators import (
#     AnnDataCatManager,
#     CatManager,
#     CurateLookup,
#     Curator,
#     DataFrameCatManager,
#     MuDataCatManager,
#     TiledbsomaCatManager,
# )
# from lamindb.models import (
#     BasicRecord,
#     CanCurate,
#     FeatureValue,
#     HasParents,
#     IsVersioned,
#     ParamValue,
#     Record,
#     Registry,
#     TracksRun,
#     TracksUpdates,
#     ValidateFields,
# )
from .. import errors as exceptions
from . import datasets, fields, loaders, subsettings, types

# from ._context import Context
from ._mapped_collection import MappedCollection
from ._settings import Settings
