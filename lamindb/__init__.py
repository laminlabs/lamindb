"""Lamin: data management for computational biologists."""

from . import _version

__version__ = _version.get_versions()["version"]
from ._logging import logger  # noqa

from . import storage  # noqa
from .database._notion import Dataset  # noqa
from ._configure import Configure  # noqa
