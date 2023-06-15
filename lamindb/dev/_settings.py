from typing import Literal

from lamin_logger import logger


class Settings:
    """Settings.

    Directly use instance `lamindb.settings` rather instantiating this class
    yourself.
    """

    def __init__(self):
        self._verbosity: int = 2  # info-level logging

    if_file_hash_exists: Literal[
        "warn_return_existing", "error", "warn_create_new"
    ] = "warn_return_existing"
    """Behavior if file hash exists (default: "warn_return_existing").

    One of `["warn_return_existing", "error", "warn_create_new"]`.

    FAQ: :doc:`/faq/ingest-same-file-twice`.
    """
    track_run_inputs_upon_load: bool = False
    """Upon load, add loaded files as the input of the current notebook run.

    FAQ: :doc:`/faq/track-runin`.
    """

    @property
    def verbosity(self) -> int:
        """Verbosity (default 2).

        - 0: only show 'error' messages
        - 1: also show 'warning' messages
        - 2: also show 'info' messages
        - 3: also show 'hint' messages
        - 4: also show detailed 'debug' messages

        This is based on Scanpy's and Django's verbosity setting.
        """
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity: int):
        self._verbosity = verbosity
        logger.set_verbosity(verbosity)


settings = Settings()
