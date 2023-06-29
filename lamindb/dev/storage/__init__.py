"""Storage.

.. autosummary::
   :toctree: .

   AnnDataAccessor
   UPath

"""
from lamindb_setup.dev.upath import UPath
from lamindb_setup.dev.upath import infer_filesystem as _infer_filesystem

from ._anndata_sizes import size_adata

try:
    from ._backed_access import AnnDataAccessor
except ImportError:
    pass
from ._object import infer_suffix, write_to_file
from .file import delete_storage, load_to_memory, store_object
