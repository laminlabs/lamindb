"""lamindb: Manage files & data.

Import the package::

   import lamindb as db

Browse the API:

.. autosummary::
   :toctree: .

   do
   model
   track
   admin
   dev
   settings
"""

__version__ = "0.3.dev1"
from . import admin  # noqa
from . import dev  # noqa
from . import do  # noqa
from . import track  # noqa
from ._settings import settings  # noqa
