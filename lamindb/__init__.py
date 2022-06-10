"""lamindb: Manage files & data.

Import the package::

   import lamindb as lndb

Query & load data as `DataFrame`, `AnnData` or `MuData`:

.. autosummary::
   :toctree: .

   query
   load

Modify data:

.. autosummary::
   :toctree: .

   ingest
   update
   delete

Inspect the schema:

.. autosummary::
   :toctree: .

   schema

Track the state of the DB:

.. autosummary::
   :toctree: .

   track.dataflow
   track.access
   track.integrity
   track.unused

Settings:

.. autosummary::
   :toctree: .

   settings

Admin tasks: Setup the database & update the schema.

.. autosummary::
   :toctree: .

   _admin
"""

__version__ = "0.3.dev1"
from . import _admin  # noqa
from . import track  # noqa
from ._delete import delete  # noqa
from ._ingest import ingest  # noqa
from ._load import load  # noqa
from ._query import query  # noqa
from ._schema import schema  # noqa
from ._settings import settings  # noqa
from ._update import update  # noqa
