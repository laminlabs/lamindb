from pathlib import Path
from typing import Literal, Mapping, Tuple, Union

import lamindb_setup as ln_setup
from lamin_utils import logger
from upath import UPath


class Settings:
    """Settings.

    Directly use instance `lamindb.settings` rather than instantiating this
    class yourself.
    """

    def __init__(self):
        self._verbosity: int = 4  # hint-level logging
        logger.set_verbosity(self._verbosity)

    upon_file_create_if_hash_exists: Literal[
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
        return ln_setup.settings.storage.root

    @storage.setter
    def storage(
        self, path_kwargs: Union[str, Path, UPath, Tuple[Union[str, UPath], Mapping]]
    ):
        if isinstance(path_kwargs, tuple):
            path, kwargs = path_kwargs
        else:
            path, kwargs = path_kwargs, {}
        ln_setup.set.storage(path, **kwargs)

    @property
    def verbosity(self) -> int:
        """Verbosity (default 4 / 'hint').

        - 0: ❌ only show 'error' messages
        - 1: ❗ also show 'warning' messages
        - 2: ✅ also show 'success' and 'save' messages
        - 3: 💡 also show 'info' messages
        - 4: 💡 also show 'hint' messages
        - 5: 🐛 also show detailed 'debug' messages

        This is based on Scanpy's and Django's verbosity setting.
        """
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity: int):
        self._verbosity = verbosity
        logger.set_verbosity(verbosity)


settings = Settings()
