from __future__ import annotations

import os
from typing import TYPE_CHECKING, Literal, Mapping

import lamindb_setup as ln_setup
from lamin_utils import logger
from lamindb_setup._set_managed_storage import set_managed_storage
from lamindb_setup.core._settings import settings as setup_settings
from lamindb_setup.core._settings_instance import sanitize_git_repo_url

from .subsettings._creation_settings import CreationSettings, creation_settings
from .subsettings._transform_settings import TransformSettings, transform_settings

if TYPE_CHECKING:
    from pathlib import Path

    from upath import UPath

VERBOSITY_TO_INT = {
    "error": 0,  # 40
    "warning": 1,  # 30
    "success": 2,  # 25
    "info": 3,  # 20
    "hint": 4,  # 15
    "debug": 5,  # 10
}
VERBOSITY_TO_STR: dict[int, str] = dict(
    [reversed(i) for i in VERBOSITY_TO_INT.items()]  # type: ignore
)


class Settings:
    """Settings.

    Use ``lamindb.settings`` instead of instantiating this class yourself.
    """

    def __init__(self, git_repo: str | None):
        self._verbosity_int: int = 1  # warning-level logging
        logger.set_verbosity(self._verbosity_int)
        self._sync_git_repo: str | None = git_repo

    @property
    def creation(self) -> CreationSettings:
        """Record creation settings.

        For example, `ln.settings.creation.search_names = False` will disable
        searching for records with similar names during creation.
        """
        return creation_settings

    track_run_inputs: bool = True
    """Track files as input upon `.load()`, `.cache()` and `.backed()`.

    Requires a global run context with :func:`~lamindb.track` was created!

    FAQ: :doc:`/faq/track-run-inputs`
    """
    __using_key: str | None = None
    _using_storage: str | None = None

    @property
    def _using_key(self) -> str | None:
        """Key for Django database settings."""
        return self.__using_key

    @_using_key.setter
    def _using_key(self, value: str | None):
        ln_setup.settings._using_key = value
        self.__using_key = value

    @property
    def _storage_settings(self) -> ln_setup.core.StorageSettings:
        if self._using_storage is None:
            storage_settings = ln_setup.settings.storage
        else:
            storage_settings = ln_setup.core.StorageSettings(root=self._using_storage)
        return storage_settings

    @property
    def transform(self) -> TransformSettings:
        """Transform settings."""
        return transform_settings

    @property
    def sync_git_repo(self) -> str | None:
        """Sync transforms with scripts in git repository.

        Provide the full git repo URL.
        """
        return self._sync_git_repo

    @sync_git_repo.setter
    def sync_git_repo(self, value) -> None:
        """Sync transforms with scripts in git repository.

        Provide the full git repo URL.
        """
        self._sync_git_repo = sanitize_git_repo_url(value)
        assert self._sync_git_repo.startswith("https://")

    @property
    def storage(self) -> Path | UPath:
        """Default storage location (a path to its root).

        Examples:

        You can switch to another managed storage location via:

        >>> ln.settings.storage = "s3://some-bucket"

        You can also pass additional fsspec kwargs via:

        >>> kwargs = dict(
        >>>     profile="some_profile", # fsspec arg
        >>>     cache_regions=True # fsspec arg for s3
        >>> )
        >>> ln.settings.storage = "s3://some-bucket", kwargs
        """
        return self._storage_settings.root

    @storage.setter
    def storage(self, path_kwargs: str | Path | UPath | tuple[str | UPath, Mapping]):
        if isinstance(path_kwargs, tuple):
            path, kwargs = path_kwargs
        else:
            path, kwargs = path_kwargs, {}
        set_managed_storage(path, **kwargs)

    @property
    def storage_local(self) -> Path:
        """An additional local default storage (a path to its root).

        Is only available if :attr:`~lamindb.setup.core.InstanceSettings.keep_artifacts_local` is enabled.

        Guide: :doc:`faq/keep-artifacts-local`

        Shortcut for: `ln.setup.settings.instance.storage_local.root`
        """
        return ln_setup.settings.instance.storage_local.root

    @storage_local.setter
    def storage_local(self, local_root: Path):
        ln_setup.settings.instance.storage_local = local_root

    @property
    def verbosity(self) -> str:
        """Logger verbosity (default 'warning').

        - 'error': ❌ only show error messages
        - 'warning': ❗ also show warning messages
        - 'success': ✅ also show success and save messages
        - 'info': 💡 also show info messages
        - 'hint': 💡 also show hint messages
        - 'debug': 🐛 also show detailed debug messages
        """
        return VERBOSITY_TO_STR[self._verbosity_int]

    @verbosity.setter
    def verbosity(self, verbosity: str | int):
        if isinstance(verbosity, str):
            verbosity_int = VERBOSITY_TO_INT[verbosity]
        else:
            verbosity_int = verbosity
        self._verbosity_int = verbosity_int
        logger.set_verbosity(verbosity_int)


if os.environ.get("LAMINDB_MULTI_INSTANCE") == "true":
    git_repo = None
else:
    git_repo = setup_settings.instance.git_repo

settings = Settings(git_repo=git_repo)
