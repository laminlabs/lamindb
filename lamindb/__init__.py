"""LaminDB: Manage data & analyses.

Import the package::

   import lamindb as ln

Browse the API:

.. autosummary::
   :toctree: .

   db
   schema
   settings
   datasets
   nb
   dev

"""

__version__ = "0.9.2"

import warnings

from lndb_setup import settings  # noqa
from lndb_setup._migrate import check_migrate as _check_migrate

from . import _check_versions  # executes checks during import

if settings.instance.storage_root is None:
    raise RuntimeError("Please run `lndb init` to configure an instance.")
_check_migrate(usettings=settings.user, isettings=settings.instance)

from . import datasets  # noqa
from . import db  # noqa
from . import dev  # noqa
from . import schema  # noqa
from ._nb import nb  # noqa

settings.__doc__ = """Settings.

This re-exports `lndb_setup.settings <https://lamin.ai/docs/lndb-setup/lndb_setup.settings>`__.
"""
