"""Storage API.

Valid suffixes.

.. autodata:: VALID_SUFFIXES

Array accessors.

.. autoclass:: AnnDataAccessor
.. autoclass:: SpatialDataAccessor
.. autoclass:: BackedAccessor
"""

from lamindb_setup.core.upath import LocalPathClasses, UPath, infer_filesystem

from ._backed_access import AnnDataAccessor, BackedAccessor, SpatialDataAccessor
from ._tiledbsoma import save_tiledbsoma_experiment
from ._valid_suffixes import VALID_SUFFIXES
from .objects import infer_suffix, write_to_disk
from .paths import delete_storage
