"""Manage data.

View data:

.. autosummary::
   :toctree: .

   view

Select & load:

.. autosummary::
   :toctree: .

   select
   load

Ingest data:

.. autosummary::
   :toctree: .

   Ingest

Modify metadata:

.. autosummary::
   :toctree: .

   insert
   update

Delete data:

.. autosummary::
   :toctree: .

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

from lndb_hub import hub  # noqa, currently not documented as being overhauled

from ..dev.db import session  # noqa
from ..dev.db._insert import insert  # noqa
from ..dev.db._link import link  # noqa
from ..dev.db._select import select  # noqa
from ..dev.db._update import update  # noqa
from ._delete import delete  # noqa
from ._ingest import Ingest  # noqa
from ._load import load  # noqa
from ._view import view  # noqa
