from typing import Literal


class CreationSettings:
    artifact_if_hash_exists: Literal[
        "warn_return_existing", "error", "warn_create_new"
    ] = "warn_return_existing"
    """Behavior if file hash exists (default `"warn_return_existing"`).

    One of `["warn_return_existing", "error", "warn_create_new"]`.

    FAQ: :doc:`/faq/idempotency`
    """
    artifact_skip_size_hash: bool = False
    """To speed up registering high numbers of files (default `False`).

    This bypasses queries for size and hash to AWS & GCP.

    It speeds up file creation by about a factor 100.
    """
    search_names: bool = True
    """To speed up creating records (default `True`).

    If `True`, search for alternative names.

    FAQ: :doc:`/faq/idempotency`
    """
    artifact_silence_missing_run_warning: bool = False
    """Silence warning about missing run & transform during artifact creation."""
    _artifact_use_virtual_keys: bool = True
    """Treat `key` parameter in :class:`~lamindb.Artifact` as virtual.

    If `True`, the `key` is **not** used to construct file paths, but file paths are
    based on the `uid` of artifact.
    """


creation_settings = CreationSettings()
