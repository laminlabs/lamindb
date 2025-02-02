"""Exceptions.

.. autosummary::
   :toctree: .

   ValidationError
   InvalidArgument
   DoesNotExist
   NotebookNotSaved
   MissingContextUID
   UpdateContext
   IntegrityError
   RecordNameChangeIntegrityError

"""

# inheriting from SystemExit has the sole purpose of suppressing
# the traceback - this isn't optimal but the current best solution
# https://laminlabs.slack.com/archives/C04A0RMA0SC/p1726856875597489


class ValidationError(SystemExit):
    """Validation error."""

    pass


class InvalidArgument(SystemExit):
    """Invalid method or function argument."""

    pass


class TrackNotCalled(SystemExit):
    """`ln.track()` wasn't called."""

    pass


class NotebookNotSaved(SystemExit):
    """Notebook wasn't saved."""

    pass


# equivalent to Django's DoesNotExist
# and SQLAlchemy's NoResultFound
class DoesNotExist(SystemExit):
    """No record found."""

    pass


class InconsistentKey(Exception):
    """Inconsistent transform or artifact `key`."""

    pass


class RecordNameChangeIntegrityError(SystemExit):
    """Custom exception for name change errors."""

    pass


class FieldValidationError(SystemExit):
    """Field validation error."""

    pass


# -------------------------------------------------------------------------------------
# run context
# -------------------------------------------------------------------------------------


class IntegrityError(Exception):
    """Integrity error.

    For instance, it's not allowed to delete artifacts outside managed storage
    locations.
    """

    pass


class MissingContextUID(SystemExit):
    """User didn't define transform settings."""

    pass


class UpdateContext(SystemExit):
    """Transform settings require update."""

    pass
