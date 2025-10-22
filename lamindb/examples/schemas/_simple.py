from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ... import Schema


def valid_features() -> Schema:
    """A `DataFrame` schema that validates that columns map on existing features.

    .. literalinclude:: scripts/define_valid_features.py
        :language: python
    """
    from ... import Schema

    try:
        return Schema.get(name="valid_features")
    except Schema.DoesNotExist:
        try:
            from . import define_valid_features  # noqa

            return Schema.get(name="valid_features")
        except Schema.DoesNotExist:
            importlib.reload(define_valid_features)
            return Schema.get(name="valid_features")
