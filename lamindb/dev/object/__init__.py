"""Object → storage.

Memory → file:

.. autosummary::
   :toctree: .

   anndata_to_h5ad
"""

from ._anndata import anndata_to_h5ad
from ._core import infer_file_suffix, write_to_file
