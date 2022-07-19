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
from nbproject import __version__ as nbproject_version
from packaging import version

if version.parse(nbproject_version) < version.parse("0.4.3"):
    raise RuntimeError("lamindb needs nbproject >= 0.4.3")
del version, nbproject_version

from ._delete import delete  # noqa
from ._ingest import ingest  # noqa
from ._load import load  # noqa
from ._query import query  # noqa
from ._update import update  # noqa
