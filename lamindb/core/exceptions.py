"""Exceptions.

The registry base class:

.. autosummary::
   :toctree: .

   ValidationError
   NotebookNotSavedError
   NoTitleError
   MissingTransformSettings
   UpdateTransformSettings

"""


class ValidationError(SystemExit):
    """Validation error: not mapped in registry."""

    pass


# -------------------------------------------------------------------------------------
# ln.track() AKA run_context
# -------------------------------------------------------------------------------------


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
