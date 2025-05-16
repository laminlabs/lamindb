from __future__ import annotations

from lamindb_setup.core.upath import VALID_COMPOSITE_SUFFIXES, VALID_SIMPLE_SUFFIXES

# add new composite suffixes like so
VALID_COMPOSITE_SUFFIXES.update(
    {
        ".vitessce.json",
        "spatialdata.zarr",
        ".ome.zarr",
    }
)
# can do the same for simple valid suffixes


class VALID_SUFFIXES:
    """Valid suffixes."""

    SIMPLE: set[str] = VALID_SIMPLE_SUFFIXES
    """Simple suffixes."""
    COMPOSITE: set[str] = VALID_COMPOSITE_SUFFIXES
    """Composite suffixes."""
