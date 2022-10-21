"""Manage data.

View data:

.. autosummary::
   :toctree: .

   view

Get & select metadata:

.. autosummary::
   :toctree: .

   select
   get

Ingest & load data objects and datasets:

.. autosummary::
   :toctree: .

   Ingest
   load

Delete data:

.. autosummary::
   :toctree: .

   add
   delete

Link metadata:

.. autosummary::
   :toctree: .

   link

Get a `SQLModel <https://sqlmodel.tiangolo.com>`__ session:

.. autosummary::
   :toctree: .

   session

"""

from ..dev.db import session  # noqa
from ..dev.db._add import add  # noqa
from ..dev.db._get import get  # noqa
from ..dev.db._link import link  # noqa
from ..dev.db._select import select  # noqa
from ._delete import delete  # noqa
from ._ingest import Ingest  # noqa
from ._load import load  # noqa
from ._view import view  # noqa
