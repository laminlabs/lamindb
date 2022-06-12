"""Administrate the database.

.. autosummary::
   :toctree: .

   get_engine
   get_database_file
   setup
   insert
   insert_if_not_exists
"""

from ._engine import get_database_file, get_engine  # noqa
from ._insert import insert, insert_if_not_exists  # noqa
from ._setup import setup  # noqa
