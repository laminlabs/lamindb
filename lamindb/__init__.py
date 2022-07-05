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
   setup
   dev
   admin
"""

__version__ = "0.0.8"
from . import admin  # noqa
from . import dev  # noqa
from . import do  # noqa
from . import model  # noqa
from . import setup  # noqa
from . import track  # noqa
