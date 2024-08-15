"""Exceptions.

The registry base class:

.. autosummary::
   :toctree: .

   DoesNotExist
   ValidationError
   NotebookNotSavedError
   NoTitleError
   MissingContext
   UpdateContext
   IntegrityError

"""


class ValidationError(SystemExit):
    """Validation error: not mapped in registry."""

    pass


# inspired by Django's DoesNotExist
# equivalent to SQLAlchemy's NoResultFound
class DoesNotExist(Exception):
    """No record found."""

    pass


# -------------------------------------------------------------------------------------
# ln.context.track() AKA context
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


class MissingContext(SystemExit):
    """User didn't define transform settings."""

    pass


class UpdateContext(SystemExit):
    """Transform settings require update."""

    pass
