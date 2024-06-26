from __future__ import annotations


class TransformSettings:
    stem_uid: str | None = None
    """Defines the version family of the transform.

    For example, all notebooks of the same family have a uid that starts with
    `"FPnfDtJz8qbE"`.

    The full uids of the notebooks in this family are of form
    `"{stem_uid}{suffix_uid}"` where the `suffix_uid` encodes the semantic
    `version`.
    """
    version: str | None = None
    """The version."""
    name: str | None = None
    """A name like a notebook or script title."""


transform_settings = TransformSettings()
