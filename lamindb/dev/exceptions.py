"""Exceptions.

The registry base class:

.. autosummary::
   :toctree: .

   ValidationError

"""


class ValidationError(Exception):
    """Validation error: not mapped in registry."""

    pass
