"""lamindb: Manage data & analyses.

Import the package::

   import lamindb as db  # data scientists working only with lamindb
   import lamindb as lndb  # software engineers working with several databases

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

__version__ = "0.3a1"
from . import admin  # noqa
from . import dev  # noqa
from . import do  # noqa
from . import model  # noqa
from . import track  # noqa
from ._settings import settings  # noqa
