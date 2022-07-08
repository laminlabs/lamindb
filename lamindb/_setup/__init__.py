"""Setup API.

A call of `setup` configures settings and logs in the user.
If there isn't one yet, it creates a database.

.. autofunction:: setup_from_cli
"""

from ._settings import Settings  # noqa
from ._settings import load_settings  # noqa
from ._setup import setup_from_cli  # noqa
