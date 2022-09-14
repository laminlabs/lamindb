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
from . import _check_versions

__version__ = "0.3.8"
from lndb_setup import settings  # noqa

from . import datasets  # noqa
from . import db  # noqa
from . import dev  # noqa
from . import schema  # noqa
from ._nb import nb  # noqa

settings.__doc__ = """Settings.

This re-exports `lndb_setup.settings <https://lamin.ai/docs/lndb-setup/lndb_setup.settings>`__.
"""
