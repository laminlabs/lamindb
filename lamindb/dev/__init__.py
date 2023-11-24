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

Functionality of data registries:

.. autosummary::
   :toctree: .

   Data
   FeatureManager
   LabelManager
   IsTree

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
   MappedDataset
"""

from lamin_utils._inspect import InspectResult
from lnschema_core.models import CanValidate, Data, HasParents, IsTree, Registry

from lamindb._query_manager import QueryManager
from lamindb._query_set import QuerySet
from lamindb.dev._feature_manager import FeatureManager
from lamindb.dev._label_manager import LabelManager

from . import _data, datasets, exceptions, fields, types  # noqa
from ._mapped_dataset import MappedDataset
from ._run_context import run_context
from ._settings import Settings
