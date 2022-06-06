"""lamindb: Manage files & data."""

from . import _version

__version__ = _version.get_versions()["version"]
from . import storage  # noqa
from ._configure import Configure  # noqa
from ._db._notion import Dataset  # noqa
from ._db._sqlite import DB  # noqa
from ._ingest import ingest  # noqa
from ._logging import logger  # noqa
