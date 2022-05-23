"""Lamin: data management for computational biologists."""

from . import _version

__version__ = _version.get_versions()["version"]
from . import storage  # noqa
from ._configure import Configure  # noqa
from ._logging import logger  # noqa
from .database._notion import Dataset  # noqa
