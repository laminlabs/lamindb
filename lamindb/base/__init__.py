"""Base library.

Is available also when no instance is setup.

Modules
-------

.. autosummary::
   :toctree: .

   uids
   types
   fields
   utils

"""

# we do not document dtypes right now because it's use is mostly internal
# dtypes are documented on the lamindb.feature document and we do not
# want to introduce another documentation page for it
from . import dtypes, fields, types, uids, utils
from .utils import deprecated, doc_args

__all__ = ["dtypes", "fields", "types", "uids", "utils"]
