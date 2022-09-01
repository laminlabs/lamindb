"""LaminDB: Manage data & analyses.

Import the package::

   import lamindb as ln

Browse the API:

.. autosummary::
   :toctree: .

   db
   schema
   datasets
   nb
   dev

To retrieve settings, use `lamindb.settings <https://lamin.ai/docs/lndb-setup/lndb_setup.settings>`__.
"""
from . import _check_versions

__version__ = "0.3.4"
from lndb_setup import settings  # noqa

from . import datasets  # noqa
from . import db  # noqa
from . import dev  # noqa
from . import schema  # noqa
from ._nb import nb  # noqa
