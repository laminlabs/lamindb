"""Storage → object.

File → memory:

.. autosummary::

   h5ad_to_anndata

Utilities:

.. autosummary::

   filepath
   local
   local_filepath
"""

from ._file import filepath, local, local_filepath
from ._h5ad import h5ad_to_anndata
