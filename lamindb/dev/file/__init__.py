"""Storage → object.

File → memory:

.. autosummary::
   :toctree: .

   h5ad_to_anndata

Store files:

.. autosummary::
   :toctree: .

   store_file
   store_png
   delete_file
"""

from ._file import delete_file, load_to_memory, store_file
from ._h5ad import h5ad_to_anndata
from ._images import store_png
