"""Storage API.

Valid suffixes.

.. autodata:: VALID_SUFFIXES

Array accessors.

.. autoclass:: AnnDataAccessor
.. autoclass:: SpatialDataAccessor
.. autoclass:: BackedAccessor
"""

from typing import TYPE_CHECKING

from lamindb_setup.core.upath import LocalPathClasses, UPath, infer_filesystem

from ._valid_suffixes import VALID_SUFFIXES

if TYPE_CHECKING:
    from ._backed_access import (  # noqa: TC004
        AnnDataAccessor,
        BackedAccessor,
        SpatialDataAccessor,
    )
    from ._tiledbsoma import save_tiledbsoma_experiment  # noqa: TC004

from .objects import infer_suffix, write_to_disk
from .paths import delete_storage

__all__ = [
    "AnnDataAccessor",
    "BackedAccessor",
    "SpatialDataAccessor",
    "LocalPathClasses",
    "UPath",
    "infer_filesystem",
    "VALID_SUFFIXES",
    "infer_suffix",
    "write_to_disk",
    "delete_storage",
    "save_tiledbsoma_experiment",
]


def __getattr__(name: str):
    """Lazy-import heavy symbols to avoid loading anndata/pyarrow at package import."""
    if (
        name == "AnnDataAccessor"
        or name == "BackedAccessor"
        or name == "SpatialDataAccessor"
    ):
        from ._backed_access import AnnDataAccessor, BackedAccessor, SpatialDataAccessor

        return (
            AnnDataAccessor
            if name == "AnnDataAccessor"
            else BackedAccessor
            if name == "BackedAccessor"
            else SpatialDataAccessor
        )
    if name == "save_tiledbsoma_experiment":
        from ._tiledbsoma import save_tiledbsoma_experiment

        return save_tiledbsoma_experiment
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
