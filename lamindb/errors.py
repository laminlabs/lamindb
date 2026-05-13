"""Errors.

Django.

.. autoexception:: ObjectDoesNotExist
.. autoexception:: MultipleObjectsReturned

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
.. autoexception:: FileNotInDevDir
.. autoexception:: BranchAlreadyExists

"""

# -------------------------------------------------------------------------------------
# Django
# -------------------------------------------------------------------------------------

from django.core.exceptions import (
    MultipleObjectsReturned,  # noqa: F401
    ObjectDoesNotExist,  # noqa: F401
)

ObjectDoesNotExist.__doc__ = """Object does not exist.

This is an alias for `django.core.exceptions.ObjectDoesNotExist`.
"""
DoesNotExist = ObjectDoesNotExist  # backward compat

MultipleObjectsReturned.__doc__ = """Multiple objects returned.

This is an alias for `django.core.exceptions.MultipleObjectsReturned`.
"""
MultipleResultsFound = MultipleObjectsReturned  # backward compat

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


class FileNotInDevDir(Exception):
    """File path is not within the configured dev directory."""

    pass


class BranchAlreadyExists(Exception):
    """Branch already exists.

    Raised when creating a branch with `ln.setup.switch(..., create=True)` and
    a branch with the given name or uid already exists. Consistent with `git switch -c`.
    """

    pass
