"""Setup API.

A call of `setup` configures settings and logs in the user.
If there isn't one yet, it creates a database.

.. autofunction:: setup

To retrieve all set and inferred settings after setup, use:

.. autofunction:: settings

..
   autosummary does not work with two objects that merely differ by capitalization

.. autoclass:: Settings
   :members:

"""

from ._settings import Settings, settings  # noqa
from ._setup import setup
