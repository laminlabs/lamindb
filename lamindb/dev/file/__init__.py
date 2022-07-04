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
"""

from ._file import store_file
from ._h5ad import h5ad_to_anndata
from ._images import store_png
