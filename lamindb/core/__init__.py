"""Core library.

Registries:

.. autosummary::
   :toctree: .

   Record
   BasicRecord
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

   BaseCurator
   DataFrameCurator
   AnnDataCurator
   MuDataCurator
   SOMACurator
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
   exceptions
   subsettings
   logger

"""

from lamin_utils import logger
from lamin_utils._inspect import InspectResult

from lamindb._query_manager import QueryManager
from lamindb._query_set import QuerySet, RecordList
from lamindb.core._feature_manager import FeatureManager, ParamManager
from lamindb.core._label_manager import LabelManager
from lamindb.curators import (
    AnnDataCurator,
    BaseCurator,
    CurateLookup,
    DataFrameCurator,
    MuDataCurator,
    SOMACurator,
)
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

from . import _data, datasets, exceptions, fields, loaders, subsettings, types
from ._context import Context
from ._mapped_collection import MappedCollection
from ._settings import Settings
