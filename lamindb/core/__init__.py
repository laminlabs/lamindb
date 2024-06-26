"""Core library.

Registries:

.. autosummary::
   :toctree: .

   Registry
   QuerySet
   QueryManager
   RecordsList
   HasFeatures
   HasParams
   FeatureManager
   ParamManager
   LabelManager
   IsVersioned
   CanValidate
   HasParents
   TracksRun
   TracksUpdates
   InspectResult
   fields

Annotators:

.. autosummary::
   :toctree: .

   DataFrameAnnotator
   AnnDataAnnotator
   MuDataAnnotator
   AnnotateLookup

Other:

.. autosummary::
   :toctree: .

   Settings
   MappedCollection
   run_context

Modules:

.. autosummary::
   :toctree: .

   datasets
   storage
   types
   exceptions
   subsettings

"""

from lamin_utils._inspect import InspectResult
from lnschema_core.models import (
    CanValidate,
    HasFeatures,
    HasParams,
    HasParents,
    IsVersioned,
    Registry,
    TracksRun,
    TracksUpdates,
)

from lamindb._annotate import (
    AnnDataAnnotator,
    AnnotateLookup,
    DataFrameAnnotator,
    MuDataAnnotator,
)
from lamindb._query_manager import QueryManager
from lamindb._query_set import QuerySet, RecordsList
from lamindb.core._feature_manager import FeatureManager, ParamManager
from lamindb.core._label_manager import LabelManager

from . import _data, datasets, exceptions, fields, subsettings, types
from ._mapped_collection import MappedCollection
from ._run_context import run_context
from ._settings import Settings
