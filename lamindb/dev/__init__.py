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
   datasets
   hashing
   storage
   Settings
   run_context
"""

from lnschema_core.models import Data, Registry, SynonymsAware, ValidationAware

from lamindb._query_manager import QueryManager
from lamindb._query_set import QuerySet
from lamindb.dev._feature_manager import FeatureManager

from .._context import run_context
from . import datasets  # noqa
from . import _data
from ._settings import Settings
