"""lamindb: Manage files & data.

Import the package::

   import lamindb as lndb

Main functionality:

.. autosummary::
   :toctree: .

   ingest
   db
   diagram

Settings:

.. autosummary::
   :toctree: .

   settings
"""

__version__ = "0.2.1"
from . import storage  # noqa
from ._db._notion import Dataset  # noqa
from ._db._sqlite import db  # noqa
from ._diagram import diagram  # noqa
from ._ingest import ingest  # noqa
from ._logging import logger  # noqa
from ._settings import settings  # noqa
