from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

from lamindb_setup.core import deprecated

if TYPE_CHECKING:
    from ... import Schema


@deprecated("Define a feature type instead.")
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


def _a_feature_type() -> Schema:
    """Load the feature-type schema used by tests."""
    from ... import Schema
    from . import define_a_feature_type

    # Inferred artifact schemas can reuse this itype, so identify the script schema by uid.
    schema = define_a_feature_type.schema
    if Schema.filter(uid=schema.uid).exists():
        return schema
    importlib.reload(define_a_feature_type)
    return define_a_feature_type.schema
