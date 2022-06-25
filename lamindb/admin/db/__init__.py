"""Administrate the database.

.. autosummary::
   :toctree: .

   get_engine
   insert
   insert_if_not_exists
"""

from ._engine import get_engine  # noqa
from ._insert import insert, insert_if_not_exists  # noqa
