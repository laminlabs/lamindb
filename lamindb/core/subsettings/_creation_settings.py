class CreationSettings:
    search_names: bool = True
    """Switch off to speed up creating records (default `True`).

    If `True`, search for alternative names and avoids duplicates.

    FAQ: :doc:`/faq/idempotency`
    """
    artifact_skip_size_hash: bool = False
    """To speed up registering high numbers of files (default `False`).

    This bypasses queries for size and hash to AWS & GCP.

    It speeds up file creation by about a factor 100.
    """
    artifact_silence_missing_run_warning: bool = False
    """Silence warning about missing run & transform during artifact creation (default `False`)."""
    _artifact_use_virtual_keys: bool = True
    """Treat `key` parameter in :class:`~lamindb.Artifact` as virtual.

    If `True`, the `key` is **not** used to construct file paths, but file paths are
    based on the `uid` of artifact.
    """


creation_settings = CreationSettings()
