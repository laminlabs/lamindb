"""Query, load and modify data.

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

Link features and metadata:

.. autosummary::
   :toctree: .

   link

Share data on the hub:

.. autosummary::
   :toctree: .

   hub
"""

from lndb_hub import hub  # noqa

from ._delete import delete  # noqa
from ._ingest import ingest  # noqa
from ._link import link  # noqa
from ._load import load  # noqa
from ._query import query  # noqa
from ._update import update  # noqa
