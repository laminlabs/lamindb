"""Base library.

Settings:

.. autosummary::
   :toctree: .

   Settings

Modules:

.. autosummary::
   :toctree: .

   types
   subsettings

Utils:

.. autosummary::
   :toctree: .

   doc_args
   deprecated

"""

from lamindb_setup.core import deprecated, doc_args

from . import subsettings, types
from ._settings import Settings
