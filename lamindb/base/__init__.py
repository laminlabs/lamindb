"""Base library.

Is available also when no instance is connected.

Modules:

.. autosummary::
   :toctree: .

   types

Utils:

.. autosummary::
   :toctree: .

   doc_args
   deprecated

"""

from lamindb_setup.core import deprecated, doc_args

from . import types
