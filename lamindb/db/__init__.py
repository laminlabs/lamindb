"""Manage data.

Query & load:

.. autosummary::
   :toctree: .

   query
   load

Ingest data:

.. autosummary::
   :toctree: .

   ingest

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

A `SQLModel <https://sqlmodel.tiangolo.com>`__ session:

.. autosummary::
   :toctree: .

   session

"""

from lndb_hub import hub  # noqa, currently not documented as being overhauled

from ..dev.db import session  # noqa
from ._delete import delete  # noqa
from ._ingest import ingest  # noqa
from ._insert import insert  # noqa
from ._link import link  # noqa
from ._load import load  # noqa
from ._query import query  # noqa
from ._update import update  # noqa
