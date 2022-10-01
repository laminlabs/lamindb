"""LaminDB: Manage data & analyses.

Import the package::

   import lamindb as ln

Browse the API:

.. autosummary::
   :toctree: .

   DB
   schema
   settings
   datasets
   nb
   dev

"""
from . import _check_versions

__version__ = "0.4.0"
from lndb_setup import settings  # noqa
from lndb_setup import _settings_store

from . import datasets  # noqa
from .db import DB  # noqa

if _settings_store.current_instance_file.exists():
    # lndb_setup has been run locally and a user wants to use
    # lamindb.db as a client that runs on this static local setting
    db = DB(settings)
else:
    # `DB` is a class that needs to be instantiated with settings
    # loaded from a user session and a selected instance
    pass

from . import dev  # noqa
from . import schema  # noqa
from ._nb import nb  # noqa

settings.__doc__ = """Settings.

This re-exports `lndb_setup.settings <https://lamin.ai/docs/lndb-setup/lndb_setup.settings>`__.
"""
del _settings_store  # clean up the auto-lookup
