"""Developer API.

The registry base class:

.. autosummary::
   :toctree: .

   Registry

Queries of registries:

.. autosummary::
   :toctree: .

   QuerySet
   QueryManager
   RecordsList

Functionality of data registries:

.. autosummary::
   :toctree: .

   Data
   FeatureManager
   LabelManager
   IsTree
   IsVersioned

Functionality of metadata registries:

.. autosummary::
   :toctree: .

   CanValidate
   HasParents
   InspectResult

Auxiliary tools:

.. autosummary::
   :toctree: .

   run_context
   datasets
   hashing
   storage
   fields
   Settings
   types
   exceptions
   MappedCollection
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
from lamindb.dev._feature_manager import FeatureManager
from lamindb.dev._label_manager import LabelManager

from . import _data, datasets, exceptions, fields, types
from ._mapped_collection import MappedCollection
from ._run_context import run_context
from ._settings import Settings
