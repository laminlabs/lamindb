"""Storage.

.. autosummary::
   :toctree: .

   AnnDataAccessor
   BackedAccessor
   UPath

"""
from lamindb_setup.dev.upath import UPath
from lamindb_setup.dev.upath import infer_filesystem as _infer_filesystem

from ._anndata_sizes import size_adata

try:
    from ._backed_access import AnnDataAccessor, BackedAccessor
except ImportError:
    pass
from .file import delete_storage, load_to_memory, store_object
from .object import infer_suffix, write_to_file
