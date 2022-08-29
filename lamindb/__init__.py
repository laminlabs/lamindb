"""lamindb: Manage data & analyses.

Import the package::

   import lamindb as db  # or lndb

Browse the API:

.. autosummary::
   :toctree: .

   do
   schema
   track
   datasets
   dev
   nb
   session

To retrieve settings, use `db.settings <https://lamin.ai/docs/lndb-setup/lndb_setup.settings>`__.
"""
from . import _check_versions

__version__ = "0.3.3"
from lndb_setup import settings  # noqa

from . import datasets  # noqa
from . import dev  # noqa
from . import do  # noqa
from . import schema  # noqa
from . import track  # noqa
from ._nb import nb  # noqa
from .dev.db import session  # noqa
