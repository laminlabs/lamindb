"""Setup API.

The settings set during setup are linked below.

In to retrieve them after setup, use the instance `setup.settings`.

.. autosummary::
    :toctree:

    Settings
    setup
    settings
"""

from ._settings import Settings, settings  # noqa
from ._setup import setup
