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

Utils
-----

.. autodecorator:: doc_args
.. autodecorator:: deprecated

"""

from lamindb_setup.core import deprecated, doc_args

from . import dtypes, fields, types, uids
