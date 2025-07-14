from __future__ import annotations

import os
import sys
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import colors, logger
from lamindb_setup import settings as setup_settings
from lamindb_setup._set_managed_storage import set_managed_storage
from lamindb_setup.core import deprecated
from lamindb_setup.core._settings_instance import sanitize_git_repo_url

from .subsettings._annotation_settings import AnnotationSettings, annotation_settings
from .subsettings._creation_settings import CreationSettings, creation_settings

if TYPE_CHECKING:
    from collections.abc import Mapping
    from pathlib import Path

    from lamindb_setup.core._settings_storage import StorageSettings
    from upath import UPath


def is_read_only_connection() -> bool:
    instance = setup_settings.instance
    if instance.dialect == "postgresql":
        db_url = instance.db
        return "read" in db_url or "public" in db_url
    return False


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

    Please use the global `ln.settings` object instead of instantiating this class yourself.
    """

    def __init__(self):
        self._verbosity_int: int = 1  # warning-level logging
        logger.set_verbosity(self._verbosity_int)
        self._sync_git_repo: str | None = None

    def __repr__(self) -> str:  # pragma: no cover
        if "sphinx" in sys.modules:
            return object.__repr__(self)

        cls_name = colors.green(self.__class__.__name__)
        verbosity_color = colors.yellow if self.verbosity == "warning" else colors.green
        verbosity_str = verbosity_color(self.verbosity)

        storage_root = self._storage_settings.root_as_str
        storage_str = colors.italic(storage_root)

        instance_str = colors.italic(self.instance_uid)
        track_color = colors.green if self.track_run_inputs else colors.yellow
        track_str = track_color(str(self.track_run_inputs))

        lines = [
            f"{cls_name}",
            f"  instance: {instance_str}",
            f"  storage: {storage_str}",
            f"  verbosity: {verbosity_str}",
            f"  track_run_inputs: {track_str}",
        ]

        if self.sync_git_repo:
            repo_name = (
                self.sync_git_repo.split("/")[-1]
                if "/" in self.sync_git_repo
                else self.sync_git_repo
            )
            lines.append(f"  sync_git_repo: {colors.italic(repo_name)}")

        return "\n".join(lines)

    @property
    def creation(self) -> CreationSettings:
        """SQLRecord creation settings.

        For example, `ln.settings.creation.search_names = False` will disable
        searching for records with similar names during creation.
        """
        return creation_settings

    @property
    def annotation(self) -> AnnotationSettings:
        """Artifact annotation settings.

        For example, `ln.settings.creation.search_names = False` will disable
        searching for records with similar names during creation.
        """
        return annotation_settings

    # note: this setting should probably be deprecated soon
    # warnings could then be filtered with a regular warning mechanism
    track_run_inputs: bool = True
    """Track run inputs (default `True`).

    If this setting is true, an artifact is recorded as run input upon `.load()`, `.cache()` & `.open()` provided :func:`~lamindb.track` was called in the current compute (Python, R) session.
    If :func:`~lamindb.track` was not called, you receive a warning message upon `.load()`, `.cache()` & `.open()`.

    If you switch this setting to `False`, you won't see the warning message anymore and no run inputs will be recorded.

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
    def sync_git_repo(self) -> str | None:
        """Sync transforms with scripts in git repository.

        Provide the full git repo URL.
        """
        if self._sync_git_repo is not None:
            return self._sync_git_repo
        elif os.environ.get("LAMINDB_MULTI_INSTANCE") == "true":
            return None
        else:
            return setup_settings.instance.git_repo

    @sync_git_repo.setter
    def sync_git_repo(self, value) -> None:
        """Sync transforms with scripts in git repository.

        For example: `ln.settings.sync_git_repo = https://github.com/laminlabs/redun-lamin`
        """
        self._sync_git_repo = sanitize_git_repo_url(value)
        if not self._sync_git_repo.startswith("https://"):  # pragma: nocover
            raise ValueError("git repository URL must start with 'https://'.")

    @property
    def storage(self) -> StorageSettings:
        """Current default storage location for writes.

        Examples:

        Retrieve the storage settings::

            ln.settings.storage
            #> StorageSettings(root='s3://my-bucket')

        Retrieve the storage root::

            ln.settings.storage.root
            #> UPath('s3://my-bucket')

        You can write artifacts to other storage locations by switching the current default storage location::

            ln.settings.storage = "s3://some-bucket"

        You can also pass additional fsspec kwargs via::

            kwargs = dict(
                profile="some_profile", # fsspec arg
                cache_regions=True # fsspec arg for s3
            )
            ln.settings.storage = "s3://some-bucket", kwargs
        """
        return self._storage_settings

    @storage.setter
    def storage(self, path_kwargs: str | Path | UPath | tuple[str | UPath, Mapping]):
        if isinstance(path_kwargs, tuple):
            path, kwargs = path_kwargs
            if isinstance(kwargs, str):
                kwargs = {"host": kwargs}
        else:
            path, kwargs = path_kwargs, {}
        set_managed_storage(path, **kwargs)

    @property
    def instance_uid(self) -> str:
        """The `uid` of the current instance."""
        return ln_setup.settings.instance.uid

    @property
    def cache_dir(self) -> UPath:
        """Cache root, a local directory to cache cloud files."""
        return ln_setup.settings.cache_dir

    @property
    def local_storage(self) -> StorageSettings:
        """An additional local default storage (a path to its root).

        Is only available if :attr:`~lamindb.setup.core.InstanceSettings.keep_artifacts_local` is enabled.

        Guide: :doc:`faq/keep-artifacts-local`
        """
        return ln_setup.settings.instance.local_storage

    @local_storage.setter
    def local_storage(self, local_root: Path):
        ln_setup.settings.instance.local_storage = local_root

    @property
    @deprecated("local_storage")
    def storage_local(self) -> StorageSettings:
        return self.local_storage

    @storage_local.setter
    @deprecated("local_storage")
    def storage_local(self, local_root_host: tuple[Path | str, str]):
        self.local_storage = local_root_host  # type: ignore

    @property
    def verbosity(self) -> str:
        """Logger verbosity (default `'warning'`).

        - `'error'`: only show error messages
        - `'warning'`: also show warning messages
        - `'success'`: also show success and save messages
        - `'info'`: also show info messages
        - `'hint'`: also show hint messages
        - `'debug'`: also show detailed debug messages
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


settings = Settings()
