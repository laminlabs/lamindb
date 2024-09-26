"""Exceptions.

.. autosummary::
   :toctree: .

   InvalidArgument
   DoesNotExist
   ValidationError
   NotebookNotSavedError
   NoTitleError
   MissingContextUID
   UpdateContext
   IntegrityError

"""

# inheriting from SystemExit has the sole purpose of suppressing
# the traceback - this isn't optimal but the current best solution
# https://laminlabs.slack.com/archives/C04A0RMA0SC/p1726856875597489


class InvalidArgument(SystemExit):
    """Invalid method or function argument."""

    pass


class TrackNotCalled(SystemExit):
    """wasn't called."""

    pass


class NotebookNotSaved(SystemExit):
    """Notebook wasn't saved."""

    pass


class ValidationError(SystemExit):
    """Validation error: not mapped in registry."""

    pass


# inspired by Django's DoesNotExist
# equivalent to SQLAlchemy's NoResultFound
class DoesNotExist(Exception):
    """No record found."""

    pass


# -------------------------------------------------------------------------------------
#  AKA context
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


class MissingContextUID(SystemExit):
    """User didn't define transform settings."""

    pass


class UpdateContext(SystemExit):
    """Transform settings require update."""

    pass
