"""lamindb: Manage data & analyses.

Import the package::

   import lamindb as db  # data scientists working only with lamindb
   import lamindb as lndb  # software engineers working with several databases

Browse the API:

.. autosummary::
   :toctree: .

   do
   schema
   track
   dev
   admin

To retrieve settings after setup via the CLI, use:

.. autofunction:: load_settings

..
   autosummary does not work with two objects that merely differ by capitalization

.. autoclass:: Settings
   :members:

"""
from nbproject import __version__ as nbproject_version
from packaging import version

if version.parse(nbproject_version) < version.parse("0.4.3"):
    raise RuntimeError("lamindb needs nbproject >= 0.4.3")
del version, nbproject_version

__version__ = "0.1.1"
from . import _setup  # noqa
from . import admin  # noqa
from . import dev  # noqa
from . import do  # noqa
from . import schema  # noqa
from . import track  # noqa
from ._setup import Settings, load_settings  # noqa
