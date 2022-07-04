"""Setup API.

A call of `setup` configures settings and logs in the user.
If there isn't one yet, it creates a database.

.. autofunction:: setup_from_cli

To retrieve settings after setup, use:

.. autofunction:: load_settings

..
   autosummary does not work with two objects that merely differ by capitalization

.. autoclass:: Settings
   :members:

"""

from ._settings import Settings  # noqa
from ._settings import load_settings  # noqa
from ._setup import setup_from_cli  # noqa
