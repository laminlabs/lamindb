"""Storage → object.

File → memory:

.. autosummary::
   :toctree: .

   h5ad_to_anndata

Store files:

.. autosummary::
   :toctree: .

   store_object
   store_png
   delete_storage
"""

__version__ = "0.4a1"  # denote a release candidate for 0.1.0 with 0.1rc1

from lamindb_setup.dev.upath import UPath
from lamindb_setup.dev.upath import infer_filesystem as _infer_filesystem

from ._file import delete_storage, load_to_memory, store_object
from ._h5ad import h5ad_to_anndata
from ._images import store_png
from ._subset import subset
from ._zarr import read_adata_zarr, write_adata_zarr
