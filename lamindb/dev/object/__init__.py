"""Object → storage.

Memory → file:

.. autosummary::
   :toctree: .

   anndata_to_h5ad
"""

from ._anndata import anndata_to_h5ad
from ._anndata_sizes import size_adata
from ._core import infer_suffix, write_to_file
from ._subset_anndata import _subset_anndata_dobject
