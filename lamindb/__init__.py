"""lamindb: Manage data & analyses.

Import the package::

   import lamindb as db  # or lndb

Browse the API:

.. autosummary::
   :toctree: .

   do
   schema
   meta
   track
   dev
   header

To retrieve settings, use `db.settings <https://lamin.ai/docs/lndb-setup/lndb_setup.settings>`__.
"""
from . import _check_versions

__version__ = "0.1.2"
from lndb_setup import settings  # noqa
from nbproject import header  # noqa

from . import dev  # noqa
from . import do  # noqa
from . import schema  # noqa
from . import track  # noqa
