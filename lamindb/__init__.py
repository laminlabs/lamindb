"""LaminDB: Manage data & analyses.

Import the package::

   import lamindb as ln
   import lamindb.schema as lns

The central data object, a wrapper for files, on-disk (`zarr`, etc.) and
in-memory objects (`DataFrame`, `AnnData`, etc.):

.. autosummary::
   :toctree: .

   DObject

Query & manipulate data:

.. autosummary::
   :toctree: .

   select
   add
   delete

Manipulate data with open session:

.. autosummary::
   :toctree: .

   Session

View DB content:

.. autosummary::
   :toctree: .

   view

Schema - entities and their relations:

.. autosummary::
   :toctree: .

   schema

Manage knowledge:

.. autosummary::
   :toctree: .

   knowledge

Track Jupyter notebooks:

.. autosummary::
   :toctree: .

   nb

Developer API:

.. autosummary::
   :toctree: .

   settings
   dev
"""

__version__ = "0.25.0"

# prints warning of python versions
from lamin_logger import py_version_warning

py_version_warning("3.8", "3.10")

from lndb_setup import settings  # noqa
from lndb_setup._migrate import check_migrate as _check_migrate

from . import _check_versions  # executes checks during import

if settings.instance.storage.root is None:
    raise RuntimeError("Please run `lndb init` to configure an instance.")
_check_migrate(usettings=settings.user, isettings=settings.instance)

from lnschema_core import DObject  # noqa

from . import dev  # noqa
from . import knowledge  # noqa
from . import schema  # noqa
from ._delete import delete  # noqa
from ._nb import nb  # noqa
from ._subset import subset
from ._view import view  # noqa
from .dev.db import Session  # noqa
from .dev.db._add import add  # noqa
from .dev.db._select import select  # noqa
from .dev.object._lazy_field import lazy

settings.__doc__ = """Settings.

This re-exports `lndb_setup.settings <https://lamin.ai/docs/lndb-setup/lndb_setup.settings>`__.
"""
