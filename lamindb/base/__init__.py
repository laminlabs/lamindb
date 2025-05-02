"""Base library.

Is available also when no instance is setup.

Modules:

.. autosummary::
   :toctree: .

   uids
   types
   fields

Utils:

.. autosummary::
   :toctree: .

   doc_args
   deprecated

"""

from lamindb_setup.core import deprecated, doc_args

from . import fields, types, uids
