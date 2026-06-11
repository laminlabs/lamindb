from __future__ import annotations

from pathlib import PurePosixPath
from typing import TYPE_CHECKING, Any, Iterable, Literal

from django.db import models
from lamin_utils import logger
from lamin_utils._base62 import decode as base62_to_int
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
    # `uid` is a concrete CharField on every subclass; declared here (annotation only,
    # so Django does not create a clashing field on this abstract base) so that the
    # `stem_uid`/`version` helpers and type-checkers can resolve its type.
    uid: str

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
        self._refresh_revises_if_stale = kwargs.pop("_refresh_revises_if_stale", False)
        # if `revises` is already a previous (non-latest) version at init, the user
        # consciously passed an old version: revise from the live head, like the inferred
        # case. This keeps the concurrent-demotion check (`revises` was the head at init
        # but got demoted before save -> raise) intact for genuinely explicit heads.
        if (
            not self._refresh_revises_if_stale
            and self._revises is not None
            and not self._revises.is_latest
        ):
            self._refresh_revises_if_stale = True
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
        new_uid = create_uid(revises=revises, version_tag=version_tag)
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


def max_version_uid_in_family(record: IsVersioned) -> str | None:
    """Return the uid with the maximum base62 version suffix in the version family.

    The max is computed in Python via base62 decoding rather than ordering by `uid`
    in the database: the version suffix is base62 (0-9 < A-Z < a-z), but locale-aware
    db collations (e.g. Postgres `en_US.UTF-8`) sort letters case-insensitively, so
    e.g. `...000Z` wrongly sorts after `...000a` and the chain can't grow past the
    `Z -> a` boundary. Considering the whole family (not just `is_latest`, which is
    per-branch and thus not globally unique, and excludes trashed records) guarantees
    a derived new uid never collides with an existing one.
    """
    return max(
        record.__class__.objects.filter(uid__startswith=record.stem_uid).values_list(
            "uid", flat=True
        ),
        key=lambda uid: base62_to_int(uid[-4:]),
        default=None,
    )


def create_uid(
    *,
    version_tag: str | None = None,
    n_full_id: int = 20,
    revises: IsVersioned | None = None,
) -> str:
    """Derive the uid for a new (possibly versioned) record.

    The new version suffix is derived from the family-wide max (computed in Python via base62
    decoding, so it's independent of db collation and never collides with higher versions on
    other branches or in the trash). `revises` is only read, never modified: choosing which
    record a new version supersedes/demotes is `save`'s job, since only it knows the branch the
    new record is created on (`is_latest` is per-branch).
    """
    if revises is not None:
        suid = revises.stem_uid
        max_uid = max_version_uid_in_family(revises) or revises.uid
        vuid = increment_base62(max_uid[-4:])
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
    return suid + vuid


def process_revises(
    revises: IsVersioned | None,
    version_tag: str | None,
    key: str | None,
    description: str | None,
    type_: type[IsVersioned],
) -> tuple[str, str, str, str]:
    if revises is not None and not isinstance(revises, type_):
        raise TypeError(f"`revises` has to be of type `{type_.__name__}`")
    uid = create_uid(
        revises=revises,
        version_tag=version_tag,
        n_full_id=type_._len_full_uid,
    )
    if revises is not None:
        if description is None:
            description = getattr(revises, "description", None)
        if key is None:
            key = revises.key
    return uid, version_tag, key, description


def _adjust_is_latest_when_deleting_is_versioned(
    objects: IsVersioned | Iterable[IsVersioned],
) -> list[int]:
    """After deleting (soft or permanent) versioned records, promote a new latest head.

    Accepts a single IsVersioned instance, a QuerySet, or a list of IsVersioned.
    is_latest is per-branch, so for every deleted head we promote the most recently
    created remaining version on its branch -- a global family max could leave that
    branch headless (or flip the wrong branch's record). Non-branch-aware models keep a
    single global head per family (branch is `None`). Returns the promoted pks.
    """
    if isinstance(objects, IsVersioned):
        objects = [objects]
    else:
        objects = list(objects)
    if not objects:
        return []
    id_list = [o.pk for o in objects]

    object_0 = objects[0]
    registry = type(object_0)
    db = getattr(object_0._state, "db", None) or "default"

    len_stem = registry._len_stem_uid

    from .sqlrecord import SQLRecord

    branch_aware = issubclass(registry, SQLRecord)
    # (family, branch) heads that were just deleted and need a successor
    groups = {
        (o.uid[:len_stem], o.branch_id if branch_aware else None)
        for o in objects
        if o.is_latest
    }
    if not groups:
        return []
    # for each deleted head, let the db pick its successor: the newest remaining version
    # in the same family on the same branch (LIMIT 1, no client-side scan of the family)
    promoted: list[tuple[int, str]] = []
    for stem, branch in groups:
        candidates = (
            registry.objects.using(db)
            .filter(uid__startswith=stem)
            .exclude(pk__in=id_list)
        )
        if branch_aware:
            # a concrete branch id inherently excludes the trash (branch_id=-1)
            candidates = candidates.filter(branch_id=branch)
        new_head = candidates.order_by("-created_at", "-id").values("pk", "uid").first()
        if new_head is not None:
            promoted.append((new_head["pk"], new_head["uid"]))
    if not promoted:
        return []
    pks = [pk for pk, _ in promoted]
    registry.objects.using(db).filter(pk__in=pks).update(is_latest=True)
    promoted_uids = [uid for _, uid in promoted]
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
