from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ... import Schema


def valid_features() -> Schema:
    """Return a schema for an AnnData with Ensembl gene IDs and valid features in obs.

    .. literalinclude:: scripts/define_valid_features.py
        :language: python
    """
    from ... import Schema

    docs_path = Path(__file__).parent.parent.parent.parent / "docs" / "scripts"
    if str(docs_path) not in sys.path:
        sys.path.append(str(docs_path))

    import define_valid_features  # noqa

    importlib.reload(define_valid_features)

    return Schema.get(name="valid_features")
