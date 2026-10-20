from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ... import Schema


def valid_features() -> Schema:
    """A `DataFrame` schema that validates that columns map on existing features.

    ::

        schema = ln.Schema(name="valid_features", itype=ln.Feature).save()
    """
    from ... import Feature, Schema

    try:
        return Schema.get(name="valid_features")
    except Schema.DoesNotExist:
        return Schema(name="valid_features", itype=Feature).save()
