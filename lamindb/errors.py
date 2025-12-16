"""Errors.

Django.

.. autoexception:: DoesNotExist
.. autoexception:: MultipleResultsFound

LaminDB.

.. autoexception:: ValidationError
.. autoexception:: InvalidArgument
.. autoexception:: NotebookNotSaved
.. autoexception:: UnknownStorageLocation
.. autoexception:: MissingContextUID
.. autoexception:: UpdateContext
.. autoexception:: IntegrityError
.. autoexception:: FieldValidationError
.. autoexception:: NoWriteAccess
.. autoexception:: BlobHashNotFound

"""

# -------------------------------------------------------------------------------------
# Django
# -------------------------------------------------------------------------------------

from django.core.exceptions import ObjectDoesNotExist as DoesNotExist  # noqa: F401

DoesNotExist.__doc__ = """Record does not exist.

This is an alias for `django.core.exceptions.ObjectDoesNotExist`.
"""


class MultipleResultsFound(Exception):
    """Multiple records found."""

    pass


# -------------------------------------------------------------------------------------
# lamindb
# -------------------------------------------------------------------------------------


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


class InconsistentKey(Exception):
    """Inconsistent transform or artifact `key`."""

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


class MissingContextUID(Exception):
    """User didn't define transform settings."""

    pass


class UpdateContext(Exception):
    """Transform settings require update."""

    pass


class BlobHashNotFound(Exception):
    """Blob hash not found in git or storage."""

    pass


# -------------------------------------------------------------------------------------
# CRUD
# -------------------------------------------------------------------------------------


class NoWriteAccess(Exception):
    """No write access to a space."""

    pass
