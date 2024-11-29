"""Core library.

Registries:

.. autosummary::
   :toctree: .

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

   BaseCurator
   DataFrameCurator
   AnnDataCurator
   MuDataCurator
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
   types
   exceptions
   subsettings
   logger

"""

from lamin_utils import logger
from lamin_utils._inspect import InspectResult
from lnschema_core.models import (
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

from lamindb._curate import (
    AnnDataCurator,
    BaseCurator,
    CurateLookup,
    DataFrameCurator,
    MuDataCurator,
)
from lamindb._query_manager import QueryManager
from lamindb._query_set import QuerySet, RecordList
from lamindb.core._feature_manager import FeatureManager, ParamManager
from lamindb.core._label_manager import LabelManager

from . import _data, datasets, exceptions, fields, loaders, subsettings, types
from ._context import Context
from ._mapped_collection import MappedCollection
from ._settings import Settings
