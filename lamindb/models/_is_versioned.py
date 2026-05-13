from __future__ import annotations

from pathlib import PurePosixPath
from typing import TYPE_CHECKING, Any, Iterable, Literal

from django.db import models
from django.db.models import Q
from lamin_utils import logger
from lamin_utils._base62 import increment_base62

from lamindb.base import uids
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

    version_tag: str | None = CharField(max_length=30, null=True, db_index=True)
    """Version tag (default `None`).

    Consider using `semantic versioning <https://semver.org>`__
    with `Python versioning <https://peps.python.org/pep-0440/>`__.
    """
    is_latest: bool = BooleanField(default=True, db_index=True)
    """Boolean flag that indicates whether a record is the latest in its version family."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        self._revises = kwargs.pop("revises", None)
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
    def version(self) -> str:
        """The version of an object.

        Defines version of an object within a family of objects characterized by the same `stem_uid`.

        Returns `.version_tag` if set, otherwise the last 4 characters of the `uid`.
        """
        return self.version_tag if self.version_tag else self.uid[-4:]  # type: ignore

    @version.setter
    def version(self, value: str | None) -> None:
        self.version_tag = value

    @property
    def versions(self) -> QuerySet:
        """Lists all records of the same version family.

        Example::

            artifact.versions.to_dataframe()       # all versions of the artifact in a dataframe
            artifact.versions.get(is_latest=True)  # the latest version of the artifact
        """
        return (
            self.__class__.connect(self._state.db)
            .filter(uid__startswith=self.stem_uid)
            .order_by("-created_at")
        )

    def _add_to_version_family(
        self, revises: IsVersioned, version_tag: str | None = None
    ):
        """Add current record to a version family.

        Args:
            revises: a record that belongs to the version family.
            version_tag: semantic version tag of the record.
        """
        old_uid = self.uid  # type: ignore
        new_uid, revises = create_uid(revises=revises, version_tag=version_tag)
        if (
            self.__class__.__name__ == "Artifact"
            and self._real_key is None
            and (self._key_is_virtual or self.key is None)
        ):
            from lamindb.core.storage.paths import auto_storage_key_from_artifact_uid

            old_path = self.path
            new_storage_key = auto_storage_key_from_artifact_uid(
                new_uid, self.suffix, self._overwrite_versions
            )
            new_path = old_path.rename(
                old_path.with_name(PurePosixPath(new_storage_key).name)
            )
            logger.success(f"updated path from {old_path} to {new_path}!")
        self.uid = new_uid
        self.version_tag = version_tag
        self.save()
        logger.success(f"updated uid from {old_uid} to {new_uid}!")


def bump_version(
    version: str,
    bump_type: str = "minor",
    behavior: Literal["prompt", "error", "ignore"] = "error",
) -> str:
    """Bumps the version number by major or minor depending on the bump_type flag.

    Args:
        version: The current version in "MAJOR" or "MAJOR.MINOR" format.
        bump_type: The type of version bump, either 'major' or 'minor'.

    Returns:
        The new version string.
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
    version_tag: str | None = None,
    n_full_id: int = 20,
    revises: IsVersioned | None = None,
) -> tuple[str, IsVersioned | None]:
    """This also updates revises in case it's not the latest version.

    This is why it returns revises.
    """
    if revises is not None:
        latest_in_family = (
            revises.__class__.objects.filter(uid__startswith=revises.stem_uid)
            .order_by("uid")
            .last()
        )
        if latest_in_family is not None and latest_in_family.uid != revises.uid:
            revises = latest_in_family
            logger.warning(
                f"didn't pass the latest version in `revises`, retrieved it: {revises}"
            )
        suid = revises.stem_uid
        vuid = increment_base62(revises.uid[-4:])  # type: ignore
    else:
        suid = uids.base62(n_full_id - 4)
        vuid = "0000"
    if version_tag is not None:
        if not isinstance(version_tag, str):
            raise ValueError(
                "`version` parameter must be `None` or `str`, e.g., '0.1', '1', '2', etc."
            )
        if revises is not None:
            if version_tag == revises.version_tag:
                raise ValueError(
                    f"Please change the version tag or leave it `None`, '{revises.version_tag}' is already taken"
                )
    return suid + vuid, revises


def process_revises(
    revises: IsVersioned | None,
    version_tag: str | None,
    key: str | None,
    description: str | None,
    type: type[IsVersioned],
) -> tuple[str, str, str, str, IsVersioned | None]:
    if revises is not None and not isinstance(revises, type):
        raise TypeError(f"`revises` has to be of type `{type.__name__}`")
    uid, revises = create_uid(
        revises=revises, version_tag=version_tag, n_full_id=type._len_full_uid
    )
    if revises is not None:
        if description is None:
            description = getattr(revises, "description", None)
        if key is None:
            key = revises.key
    return uid, version_tag, key, description, revises


def _adjust_is_latest_when_deleting_is_versioned(
    objects: IsVersioned | Iterable[IsVersioned],
) -> list[int]:
    """After deleting (soft or permanent) versioned records, promote new latest per version family.

    Accepts a single IsVersioned instance, a QuerySet, or a list of IsVersioned.
    Runs in 1 query (candidates + update) when objects are passed; no extra query for uids.
    Returns the list of pks that were promoted to is_latest (for testing).
    """
    if isinstance(objects, IsVersioned):
        objects = [objects]
    else:
        objects = list(objects)
    if not objects:
        return []
    id_list = [o.pk for o in objects]
    stem_uids = list({o.uid[: o._len_stem_uid] for o in objects if o.is_latest})
    if not stem_uids:
        return []
    registry = type(objects[0])
    db = getattr(objects[0]._state, "db", None) or "default"
    len_stem = registry._len_stem_uid
    # All candidates: same family as any stem_uid, not in trash and not about to be deleted
    q = Q()
    for s in stem_uids:
        q |= Q(uid__startswith=s)
    qs = registry.objects.using(db).filter(q).exclude(pk__in=id_list)
    from .sqlrecord import SQLRecord

    if issubclass(registry, SQLRecord):
        qs = qs.exclude(branch_id=-1)
    candidates = list(qs.values("pk", "uid", "created_at"))
    # per stem_uid, pick candidate with max created_at
    by_stem: dict[str, dict[str, Any]] = {}
    for c in candidates:
        stem = c["uid"][:len_stem]
        if stem not in by_stem or c["created_at"] > by_stem[stem]["created_at"]:
            by_stem[stem] = c
    if not by_stem:
        return []
    pks = [by_stem[s]["pk"] for s in by_stem]
    registry.objects.using(db).filter(pk__in=pks).update(is_latest=True)
    if pks:
        promoted_uids = [by_stem[s]["uid"] for s in by_stem]
        if len(promoted_uids) == 1:
            logger.important_hint(
                f"new latest {registry.__name__} version is: {promoted_uids[0]}"
            )
        else:
            logger.important_hint(
                f"new latest {registry.__name__} versions: {promoted_uids}"
            )
    return pks


def reconcile_is_latest_within_branch(
    registry: type[IsVersioned],
    *,
    branch_id: int,
    db: str = "default",
) -> int:
    """Keep a single is_latest=True per version family in a branch.

    Winner selection is based on newest created_at, tie-broken by highest pk.
    Returns the number of records demoted from is_latest=True to False.
    """
    len_stem = registry._len_stem_uid
    latest_records = list(
        registry.objects.using(db)
        .filter(branch_id=branch_id, is_latest=True)
        .values("pk", "uid", "created_at")
        .order_by("uid", "created_at", "pk")
    )
    if not latest_records:
        return 0
    winners_by_stem: dict[str, dict[str, Any]] = {}
    losers: list[int] = []
    for record in latest_records:
        stem = record["uid"][:len_stem]
        winner = winners_by_stem.get(stem)
        if winner is None:
            winners_by_stem[stem] = record
            continue
        if (record["created_at"], record["pk"]) > (winner["created_at"], winner["pk"]):
            losers.append(winner["pk"])
            winners_by_stem[stem] = record
        else:
            losers.append(record["pk"])
    if not losers:
        return 0
    return registry.objects.using(db).filter(pk__in=losers).update(is_latest=False)
