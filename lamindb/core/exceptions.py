"""Exceptions.

The registry base class:

.. autosummary::
   :toctree: .

   ValidationError

"""


class ValidationError(SystemExit):
    """Validation error: not mapped in registry."""

    pass
