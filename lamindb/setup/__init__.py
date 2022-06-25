"""Setup API.

The settings set during setup are linked below.

In to retrieve them after setup, use the instance `setup.settings`.

.. autosummary::
   :toctree:

   Settings
   setup
"""
from ._settings import Settings, _load  # noqa
from ._setup import setup  # noqa


# see this for context: https://github.com/laminlabs/nbproject/blob/47ec6646679347ff58d53d969294333749c2a245/nbproject/__init__.py#L64  # noqa
def __getattr__(name):
    if name == "settings":
        return _load()

    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
