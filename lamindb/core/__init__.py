"""Core library.

Registries:

.. autosummary::
   :toctree: .

   Registry
   QuerySet
   QueryManager
   RecordsList
   Data
   FeatureManager
   LabelManager
   IsTree
   IsVersioned
   CanValidate
   HasParents
   InspectResult
   fields

Classes:

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

"""

from lamin_utils._inspect import InspectResult
from lnschema_core.models import (
    CanValidate,
    Data,
    HasParents,
    IsTree,
    IsVersioned,
    Registry,
)

from lamindb._query_manager import QueryManager
from lamindb._query_set import QuerySet, RecordsList
from lamindb.core._feature_manager import FeatureManager
from lamindb.core._label_manager import LabelManager

from . import _data, datasets, exceptions, fields, types
from ._mapped_collection import MappedCollection
from ._run_context import run_context
from ._settings import Settings
