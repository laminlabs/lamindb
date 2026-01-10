"""Base library.

Is available also when no instance is setup.

Modules
-------

.. autosummary::
   :toctree: .

   uids
   types
   fields
   dtypes
   utils

"""

from . import dtypes, fields, ids, types, uids, utils
from .utils import deprecated, doc_args
