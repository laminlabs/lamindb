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

Share data on the hub:

.. autosummary::
   :toctree: .

   hub
"""

from ._delete import delete  # noqa
from ._ingest import ingest  # noqa
from ._load import load  # noqa
from ._push import hub  # noqa
from ._query import query  # noqa
from ._update import update  # noqa
