"""Exceptions.

The registry base class:

.. autosummary::
   :toctree: .

   ValidationError
   NotebookNotSavedError
   NoTitleError
   MissingTransformSettings
   UpdateTransformSettings
   IntegrityError

"""


class ValidationError(SystemExit):
    """Validation error: not mapped in registry."""

    pass


# -------------------------------------------------------------------------------------
# ln.track() AKA run_context
# -------------------------------------------------------------------------------------


class IntegrityError(Exception):
    """Integrity error.

    For instance, it's not allowed to delete artifacts outside managed storage
    locations.
    """

    pass


class NotebookNotSavedError(Exception):
    """Notebook wasn't saved."""

    pass


class NoTitleError(Exception):
    """Notebook has no title."""

    pass


class MissingTransformSettings(SystemExit):
    """User didn't define transform settings."""

    pass


class UpdateTransformSettings(SystemExit):
    """Transform settings require update."""

    pass
