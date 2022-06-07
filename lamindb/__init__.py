"""lamindb: Manage files & data.

Import the package::

   import lamindb as lndb

Main functionality:

.. autosummary::
   :toctree: .

   ingest
   DB
"""

__version__ = "0.1.0"
from . import storage  # noqa
from ._configure import Configure  # noqa
from ._db._notion import Dataset  # noqa
from ._db._sqlite import DB  # noqa
from ._ingest import ingest  # noqa
from ._logging import logger  # noqa
