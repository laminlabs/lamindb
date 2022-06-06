"""lamindb: Manage files & data.

Import the package::

   import lamindb as lndb

Main functionality:

.. autosummary::
   :toctree: .

   ingest
   DB
"""

from . import _version

__version__ = _version.get_versions()["version"]
from . import storage  # noqa
from ._configure import Configure  # noqa
from ._db._notion import Dataset  # noqa
from ._db._sqlite import DB  # noqa
from ._ingest import ingest  # noqa
from ._logging import logger  # noqa
