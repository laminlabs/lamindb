"""Developer API.

.. autosummary::
   :toctree: .

   Registry
   Data
   QuerySet
   QueryManager
   FeatureManager
   ValidationAware
   SynonymsAware
   ParentsAware
   InspectResult
   datasets
   hashing
   storage
   fields
   Settings
   run_context
   exc.ValidationError
"""

from lamin_utils._inspect import InspectResult
from lnschema_core.models import (
    Data,
    ParentsAware,
    Registry,
    SynonymsAware,
    ValidationAware,
)

from lamindb._query_manager import QueryManager
from lamindb._query_set import QuerySet
from lamindb.dev._feature_manager import FeatureManager

from . import datasets  # noqa
from . import _data, exc, fields
from ._run_context import run_context
from ._settings import Settings
