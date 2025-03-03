# ruff: noqa: TC004

from typing import (
    TYPE_CHECKING,  # noqa: F401
    overload,
)

from django.db import models
from lamin_utils import logger
from upath import UPath

from lamindb.base.fields import (
    BooleanField,
    CharField,
)

from ..base.types import (
    FieldAttr,
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
    def versions(self) -> "QuerySet":
        """Lists all records of the same version family.

        >>> new_artifact = ln.Artifact(df2, revises=artifact).save()
        >>> new_artifact.versions()
        """
        db = self._state.db
        if db is not None and db != "default":
            return self.__class__.using(db).filter(uid__startswith=self.stem_uid)  # type: ignore
        else:
            return self.__class__.filter(uid__startswith=self.stem_uid)  # type: ignore

    def _add_to_version_family(
        self, revises: "IsVersioned", version: str | None = None
    ):
        """Add current record to a version family.

        Args:
            revises: a record that belongs to the version family.
            version: semantic version of the record.
        """
        from ..core.versioning import create_uid, get_new_path_from_uid

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


# -------------------------------------------------------------------------------------
# A note on required fields at the Record level
#
# As Django does most of its validation on the Form-level, it doesn't offer functionality
# for validating the integrity of an Record object upon instantation (similar to pydantic)
#
# For required fields, we define them as commonly done on the SQL level together
# with a validator in Record (validate_required_fields)
#
# This goes against the Django convention, but goes with the SQLModel convention
# (Optional fields can be null on the SQL level, non-optional fields cannot)
#
# Due to Django's convention where CharFieldAttr has pre-configured (null=False, default=""), marking
# a required field necessitates passing `default=None`. Without the validator it would trigger
# an error at the SQL-level, with it, it triggers it at instantiation

# -------------------------------------------------------------------------------------
# A note on class and instance methods of core Record
#
# All of these are defined and tested within lamindb, in files starting with _{orm_name}.py

# -------------------------------------------------------------------------------------
# A note on maximal lengths of char fields
#
# 100 characters:
#     "Raindrops pitter-pattered on the windowpane, blurring the"
#     "city lights outside, curled up with a mug."
# A good maximal length for a name (title).
#
# 150 characters: We choose this for name maximal length because some users like long names.
#
# 255 characters:
#     "In creating a precise 255-character paragraph, one engages in"
#     "a dance of words, where clarity meets brevity. Every syllable counts,"
#     "illustrating the skill in compact expression, ensuring the essence of the"
#     "message shines through within the exacting limit."
# This is a good maximal length for a description field.


class LinkORM:
    pass


def deferred_attribute__repr__(self):
    return f"FieldAttr({self.field.model.__name__}.{self.field.name})"


FieldAttr.__repr__ = deferred_attribute__repr__  # type: ignore
