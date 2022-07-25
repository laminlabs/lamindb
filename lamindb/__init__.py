"""lamindb: Manage data & analyses.

Import the package::

   import lamindb as db  # or lndb

Browse the API:

.. autosummary::
   :toctree: .

   do
   schema
   meta
   track
   dev

To retrieve settings after setup via the CLI, use:

.. autofunction:: load_or_create_instance_settings
.. autofunction:: load_or_create_user_settings

..
   autosummary does not work with two objects that merely differ by capitalization

.. autoclass:: InstanceSettings
   :members:
.. autoclass:: UserSettings
   :members:

"""
from . import _check_versions

__version__ = "0.1.2"
from lndb_cli import (  # noqa
    InstanceSettings,
    UserSettings,
    load_or_create_instance_settings,
    load_or_create_user_settings,
)

from . import dev  # noqa
from . import do  # noqa
from . import schema  # noqa
from . import track  # noqa
