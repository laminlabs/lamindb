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
"""

from ._delete import delete  # noqa
from ._ingest import ingest  # noqa
from ._load import load  # noqa
from ._push import push_dobject, push_instance, unpush_instance  # noqa
from ._query import query  # noqa
from ._update import update  # noqa
