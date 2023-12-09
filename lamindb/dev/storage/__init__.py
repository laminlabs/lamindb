"""Storage tools.

.. autosummary::
   :toctree: .

   AnnDataAccessor
   BackedAccessor
"""
from lamindb_setup.dev.upath import LocalPathClasses, UPath, infer_filesystem

from ._anndata_sizes import size_adata
from ._backed_access import AnnDataAccessor, BackedAccessor
from .file import delete_storage, load_to_memory, store_object
from .object import infer_suffix, write_to_file
