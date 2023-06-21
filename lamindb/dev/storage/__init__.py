"""Storage.

.. autosummary::
   :toctree: .

   AnnDataAccessor
   UPath

"""
from lamindb_setup.dev.upath import UPath
from lamindb_setup.dev.upath import infer_filesystem as _infer_filesystem

from ._file import delete_storage, load_to_memory, store_object
from ._zarr import read_adata_zarr, write_adata_zarr
from .object._backed_access import AnnDataAccessor
