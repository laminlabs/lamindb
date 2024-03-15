import os
from pathlib import Path
from typing import Dict, Literal, Mapping, Optional, Tuple, Union

import lamindb_setup as ln_setup
from lamin_utils import logger
from lamindb_setup._add_remote_storage import switch_default_storage
from lamindb_setup.core._settings import settings as setup_settings
from lamindb_setup.core._settings_instance import sanitize_git_repo_url
from upath import UPath

from ._transform_settings import TransformSettings, transform

VERBOSITY_TO_INT = {
    "error": 0,  # 40
    "warning": 1,  # 30
    "success": 2,  # 25
    "info": 3,  # 20
    "hint": 4,  # 15
    "debug": 5,  # 10
}
VERBOSITY_TO_STR: Dict[int, str] = dict(
    [reversed(i) for i in VERBOSITY_TO_INT.items()]  # type: ignore
)


class Settings:
    """Settings.

    Use ``lamindb.settings`` instead of instantiating this class yourself.
    """

    def __init__(self, git_repo: Optional[str]):
        self._verbosity_int: int = 1  # warning-level logging
        logger.set_verbosity(self._verbosity_int)
        self._sync_git_repo: Optional[str] = git_repo

    upon_artifact_create_if_hash_exists: Literal[
        "warn_return_existing", "error", "warn_create_new"
    ] = "warn_return_existing"
    """Behavior if file hash exists (default `"warn_return_existing"`).

    One of `["warn_return_existing", "error", "warn_create_new"]`.

    FAQ: :doc:`/faq/idempotency`
    """
    upon_file_create_skip_size_hash: bool = False
    """To speed up registering high numbers of files (default `False`).

    This bypasses queries for size and hash to AWS & GCP.

    It speeds up file creation by about a factor 100.
    """
    upon_create_search_names: bool = True
    """To speed up creating Registry objects (default `True`).

    If `True`, search for alternative names.

    FAQ: :doc:`/faq/idempotency`
    """
    track_run_inputs: bool = True
    """Track files as input upon `.load()`, `.stage()` and `.backed()`.

    Requires a global run context with :func:`~lamindb.track` was created!

    FAQ: :doc:`/faq/track-run-inputs`
    """
    silence_file_run_transform_warning: bool = False
    """Silence warning about missing run & transform during file creation."""
    artifact_use_virtual_keys: bool = True
    """Treat `key` parameter in :class:`~lamindb.Artifact` as virtual.

    If `True`, the `key` is **not** used to construct file paths, but file paths are
    based on the `uid` of artifact.
    """
    __using_key: Optional[str] = None
    _using_storage: Optional[str] = None

    @property
    def _using_key(self) -> Optional[str]:
        """Key for Django database settings."""
        return self.__using_key

    @_using_key.setter
    def _using_key(self, value: Optional[str]):
        ln_setup.settings._using_key = value
        self.__using_key = value

    @property
    def _storage_settings(self) -> ln_setup.dev.StorageSettings:
        if self._using_storage is None:
            storage_settings = ln_setup.settings.storage
        else:
            storage_settings = ln_setup.dev.StorageSettings(root=self._using_storage)
        return storage_settings

    @property
    def transform(self) -> TransformSettings:
        """Transform settings."""
        return transform

    @property
    def sync_git_repo(self) -> Optional[str]:
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
    def storage(self) -> Union[Path, UPath]:
        """Default storage location (a path to its root).

        Examples:

        You can set the root via:

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
    def storage(
        self, path_kwargs: Union[str, Path, UPath, Tuple[Union[str, UPath], Mapping]]
    ):
        logger.warning(
            "you'll no longer be able to set arbitrary storage locations soon"
        )
        if isinstance(path_kwargs, tuple):
            path, kwargs = path_kwargs
        else:
            path, kwargs = path_kwargs, {}
        switch_default_storage(path, **kwargs)

    @property
    def verbosity(self) -> str:
        """Logger verbosity (default 'warning').

        - 'error': ❌ only show error messages
        - 'warning': ❗ also show warning messages
        - 'success': ✅ also show success and save messages
        - 'info': 💡 also show info messages
        - 'hint': 💡 also show hint messages
        - 'debug': 🐛 also show detailed debug messages

        This is based on Scanpy's and Django's verbosity setting.
        """
        return VERBOSITY_TO_STR[self._verbosity_int]

    @verbosity.setter
    def verbosity(self, verbosity: Union[str, int]):
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
