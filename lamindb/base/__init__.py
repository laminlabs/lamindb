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

from . import dtypes, fields, types, uids, utils
from .utils import deprecated, doc_args

__all__ = ["dtypes", "fields", "types", "uids", "utils"]
