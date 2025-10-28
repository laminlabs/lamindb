"""Errors.

.. autoexception:: ValidationError
.. autoexception:: InvalidArgument
.. autoexception:: DoesNotExist
.. autoexception:: NotebookNotSaved
.. autoexception:: UnknownStorageLocation
.. autoexception:: MissingContextUID
.. autoexception:: UpdateContext
.. autoexception:: IntegrityError
.. autoexception:: FieldValidationError
.. autoexception:: SQLRecordNameChangeIntegrityError
.. autoexception:: NoWriteAccess
.. autoexception:: BlobHashNotFound

"""


# inheriting from SystemExit has the sole purpose of suppressing
# the traceback - this isn't optimal but the current best solution
# https://laminlabs.slack.com/archives/C04A0RMA0SC/p1726856875597489


class ValidationError(Exception):
    """Validation error."""

    pass


class InvalidArgument(Exception):
    """Invalid method or function argument."""

    pass


class TrackNotCalled(Exception):
    """`ln.track()` wasn't called."""

    pass


class NotebookNotSaved(Exception):
    """Notebook wasn't saved."""

    pass


class UnknownStorageLocation(Exception):
    """Path is not contained in any known storage location."""

    pass


class NoStorageLocationForSpace(Exception):
    """No storage location found for space."""

    pass


# equivalent to Django's DoesNotExist
# and SQLAlchemy's NoResultFound
class DoesNotExist(Exception):
    """No record found."""

    pass


class MultipleResultsFound(Exception):
    """Multiple records found."""

    pass


class InconsistentKey(Exception):
    """Inconsistent transform or artifact `key`."""

    pass


class SQLRecordNameChangeIntegrityError(Exception):
    """Custom exception for name change errors."""

    pass


class FieldValidationError(Exception):
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


class BlobHashNotFound(Exception):
    """Blob hash not found in git or storage."""

    pass


# -------------------------------------------------------------------------------------
# record
# -------------------------------------------------------------------------------------


class NoWriteAccess(Exception):
    """No write access to a space."""

    pass
