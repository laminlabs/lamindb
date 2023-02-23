"""LaminDB: Manage data & analyses.

Import the package::

   import lamindb as ln
   import lamindb.schema as lns

The central data object, a wrapper for files, on-disk (`zarr`, etc.) and
in-memory objects (`DataFrame`, `AnnData`, etc.):

.. autosummary::
   :toctree: .

   DObject

Data folder, a group of data objects:

.. autosummary::
   :toctree: .

   DFolder

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

__version__ = "0.29.0"  # denote a release candidate for 0.1.0 with 0.1rc1

# prints warning of python versions
from lamin_logger import py_version_warning

py_version_warning("3.8", "3.10")

from lndb import settings as _setup_settings
from lndb._migrate import check_deploy_migration as _check_migrate

from . import _check_versions  # executes checks during import

if _setup_settings.instance.storage.root is None:
    raise RuntimeError(
        "Please run `lndb load` or `lndb init` to configure an instance."
    )
_check_migrate(usettings=_setup_settings.user, isettings=_setup_settings.instance)

from lnschema_core import DFolder  # noqa
from lnschema_core import DObject  # noqa

from . import dev  # noqa
from . import knowledge  # noqa
from . import schema  # noqa
from ._delete import delete  # noqa
from ._nb import nb  # noqa
from ._settings import settings
from ._subset import subset
from ._view import view  # noqa
from .dev.db import Session  # noqa
from .dev.db._add import add  # noqa
from .dev.db._select import select  # noqa
from .dev.object._lazy_field import lazy
