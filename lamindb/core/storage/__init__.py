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
from .paths import delete_storage

if TYPE_CHECKING:
    from ._anndata_accessor import AnnDataAccessor
    from ._backed_access import BackedAccessor
    from ._spatialdata_accessor import SpatialDataAccessor
    from ._tiledbsoma import save_tiledbsoma_experiment
    from .objects import infer_suffix, write_to_disk


__all__ = [
    "AnnDataAccessor",
    "BackedAccessor",
    "LocalPathClasses",
    "SpatialDataAccessor",
    "UPath",
    "VALID_SUFFIXES",
    "delete_storage",
    "infer_filesystem",
    "infer_suffix",
    "save_tiledbsoma_experiment",
    "write_to_disk",
]

_LAZY_EXPORTS = frozenset(
    {
        "AnnDataAccessor",
        "BackedAccessor",
        "SpatialDataAccessor",
        "infer_suffix",
        "save_tiledbsoma_experiment",
        "write_to_disk",
    }
)


def __getattr__(name: str):
    if name not in _LAZY_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    if name == "AnnDataAccessor":
        from ._anndata_accessor import AnnDataAccessor as attr
    elif name == "BackedAccessor":
        from ._backed_access import BackedAccessor as attr
    elif name == "SpatialDataAccessor":
        from ._spatialdata_accessor import SpatialDataAccessor as attr
    elif name == "save_tiledbsoma_experiment":
        from ._tiledbsoma import save_tiledbsoma_experiment as attr
    else:
        from .objects import infer_suffix, write_to_disk

        attr = infer_suffix if name == "infer_suffix" else write_to_disk

    globals()[name] = attr
    return attr
