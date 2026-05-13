"""Storage-related type definitions."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from anndata import AnnData
    from mudata import MuData
    from spatialdata import SpatialData

    ScverseDataStructures = AnnData | MuData | SpatialData
else:
    # AnnData | MuData | SpatialData; Any required for union with DataFrame in objects.py
    ScverseDataStructures = Any
