"""LaminDB: Manage data & analyses.

Import the package::

   import lamindb as ln
   import lamindb.schema as lns

Query & inspect data:

.. autosummary::
   :toctree: .

   select
   view

Manipulate data:

.. autosummary::
   :toctree: .

   record
   add
   load
   subset
   delete

Schema - entities and their relations:

.. autosummary::
   :toctree: .

   schema

Knowledge management:

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

   session
   settings
   dev

"""

__version__ = "0.18.7"

import warnings

from lndb_setup import settings  # noqa
from lndb_setup._migrate import check_migrate as _check_migrate

from . import _check_versions  # executes checks during import

if settings.instance.storage_root is None:
    raise RuntimeError("Please run `lndb init` to configure an instance.")
_check_migrate(usettings=settings.user, isettings=settings.instance)

from . import dev  # noqa
from . import knowledge  # noqa
from . import nb  # noqa
from . import schema  # noqa
from ._delete import delete  # noqa
from ._load import load  # noqa
from ._record import record  # noqa
from ._subset import subset
from ._view import view  # noqa
from .dev.db import session  # noqa
from .dev.db._add import add  # noqa
from .dev.db._select import select  # noqa
from .dev.object._lazy_field import lazy

settings.__doc__ = """Settings.

This re-exports `lndb_setup.settings <https://lamin.ai/docs/lndb-setup/lndb_setup.settings>`__.
"""
