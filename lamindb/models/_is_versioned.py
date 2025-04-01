from __future__ import annotations

from typing import TYPE_CHECKING, Literal, overload

from django.db import models
from lamin_utils import logger
from lamin_utils._base62 import increment_base62
from lamindb_setup.core.upath import LocalPathClasses, UPath

from lamindb.base import ids
from lamindb.base.fields import (
    BooleanField,
    CharField,
)

if TYPE_CHECKING:  # noqa
    from lamindb.models.query_set import QuerySet


class IsVersioned(models.Model):
    """Base class for versioned models."""

    class Meta:
        abstract = True

    _len_stem_uid: int

    version: str | None = CharField(max_length=30, null=True, db_index=True)
    """Version (default `None`).

    Defines version of a family of records characterized by the same `stem_uid`.

    Consider using `semantic versioning <https://semver.org>`__
    with `Python versioning <https://peps.python.org/pep-0440/>`__.
    """
    is_latest: bool = BooleanField(default=True, db_index=True)
    """Boolean flag that indicates whether a record is the latest in its version family."""

    @overload
    def __init__(self): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        self._revises = kwargs.pop("revises") if "revises" in kwargs else None
        super().__init__(*args, **kwargs)

    @property
    def stem_uid(self) -> str:
        """Universal id characterizing the version family.

        The full uid of a record is obtained via concatenating the stem uid and version information::

            stem_uid = random_base62(n_char)  # a random base62 sequence of length 12 (transform) or 16 (artifact, collection)
            version_uid = "0000"  # an auto-incrementing 4-digit base62 number
            uid = f"{stem_uid}{version_uid}"  # concatenate the stem_uid & version_uid

        """
        return self.uid[: self._len_stem_uid]  # type: ignore

    @property
    def versions(self) -> QuerySet:
        """Lists all records of the same version family.

        >>> new_artifact = ln.Artifact(df2, revises=artifact).save()
        >>> new_artifact.versions()
        """
        return (
            self.__class__.using(self._state.db)
            .filter(uid__startswith=self.stem_uid)
            .order_by("-created_at")
        )

    def _add_to_version_family(self, revises: IsVersioned, version: str | None = None):
        """Add current record to a version family.

        Args:
            revises: a record that belongs to the version family.
            version: semantic version of the record.
        """
        old_uid = self.uid  # type: ignore
        new_uid, revises = create_uid(revises=revises, version=version)
        if self.__class__.__name__ == "Artifact" and self._key_is_virtual:
            old_path = self.path
            new_path = get_new_path_from_uid(
                old_path=old_path, old_uid=old_uid, new_uid=new_uid
            )
            new_path = UPath(old_path).rename(new_path)
            logger.success(f"updated path from {old_path} to {new_path}!")
        self.uid = new_uid
        self.version = version
        self.save()
        logger.success(f"updated uid from {old_uid} to {new_uid}!")


def bump_version(
    version: str,
    bump_type: str = "minor",
    behavior: Literal["prompt", "error", "ignore"] = "error",
) -> str:
    """Bumps the version number by major or minor depending on the bump_type flag.

    Parameters:
    version (str): The current version in "MAJOR" or "MAJOR.MINOR" format.
    bump_type (str): The type of version bump, either 'major' or 'minor'.

    Returns:
    str: The new version string.
    """
    try:
        # Split the version into major and minor parts if possible
        parts = version.split(".")
        major = int(parts[0])
        minor = int(parts[1]) if len(parts) > 1 else 0

        if bump_type == "major":
            # Bump the major version and reset the minor version
            new_version = f"{major + 1}"
        elif bump_type == "minor":
            # Bump the minor version
            new_version = f"{major}.{minor + 1}"
        else:
            raise ValueError("bump_type must be 'major' or 'minor'")

    except (ValueError, IndexError):
        if behavior == "prompt":
            new_version = input(
                f"The current version is '{version}' - please type the new version: "
            )
        elif behavior == "error":
            raise ValueError(
                "Cannot auto-increment non-integer castable version, please provide"
                " manually"
            ) from None
        else:
            logger.warning("could not auto-increment version, fix '?' manually")
            new_version = "?"
    return new_version


def set_version(version: str | None = None, previous_version: str | None = None):
    """(Auto-) set version.

    If `version` is `None`, returns the stored version.
    Otherwise sets the version to the passed version.

    Args:
        version: Version string.
        previous_version: Previous version string.
    """
    if version is None and previous_version is not None:
        version = bump_version(previous_version, bump_type="major")
    return version


def create_uid(
    *,
    version: str | None = None,
    n_full_id: int = 20,
    revises: IsVersioned | None = None,
) -> tuple[str, IsVersioned | None]:
    """This also updates revises in case it's not the latest version.

    This is why it returns revises.
    """
    if revises is not None:
        if not revises.is_latest:
            # need one more request
            revises = revises.__class__.objects.get(
                is_latest=True, uid__startswith=revises.stem_uid
            )
            logger.warning(
                f"didn't pass the latest version in `revises`, retrieved it: {revises}"
            )
        suid = revises.stem_uid
        vuid = increment_base62(revises.uid[-4:])  # type: ignore
    else:
        suid = ids.base62(n_full_id - 4)
        vuid = "0000"
    if version is not None:
        if not isinstance(version, str):
            raise ValueError(
                "`version` parameter must be `None` or `str`, e.g., '0.1', '1', '2', etc."
            )
        if revises is not None:
            if version == revises.version:
                raise ValueError(
                    f"Please increment the previous version: '{revises.version}'"
                )
    return suid + vuid, revises


def get_new_path_from_uid(old_path: UPath, old_uid: str, new_uid: str):
    if isinstance(old_path, LocalPathClasses):
        # for local path, the rename target must be full path
        new_path = old_path.as_posix().replace(old_uid, new_uid)
    else:
        # for cloud path, the rename target must be the last part of the path
        new_path = old_path.name.replace(old_uid, new_uid)
    return new_path


def process_revises(
    revises: IsVersioned | None,
    version: str | None,
    key: str | None,
    description: str | None,
    type: type[IsVersioned],
) -> tuple[str, str, str, str, IsVersioned | None]:
    if revises is not None and not isinstance(revises, type):
        raise TypeError(f"`revises` has to be of type `{type.__name__}`")
    uid, revises = create_uid(
        revises=revises, version=version, n_full_id=type._len_full_uid
    )
    if revises is not None:
        if description is None:
            description = revises.description
        if key is None:
            key = revises.key
    return uid, version, key, description, revises
