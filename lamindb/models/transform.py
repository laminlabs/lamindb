from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, overload

from django.db import models
from django.db.models import CASCADE, PROTECT, Q
from lamin_utils import logger
from lamindb_setup.core.hashing import HASH_LENGTH, hash_file, hash_string

from lamindb.base import deprecated
from lamindb.base.fields import (
    CharField,
    DateTimeField,
    ForeignKey,
    TextField,
)
from lamindb.base.users import current_user_id

from ..models._is_versioned import process_revises
from ._is_versioned import IsVersioned, _adjust_is_latest_when_deleting_is_versioned
from .run import Run, User
from .sqlrecord import (
    BaseSQLRecord,
    IsLink,
    SQLRecord,
    init_self_from_db,
    update_attributes,
)

if TYPE_CHECKING:
    from datetime import datetime
    from pathlib import Path

    from lamindb.base.types import TransformKind

    from .artifact import Artifact
    from .block import TransformBlock
    from .project import Project, Reference
    from .query_manager import RelatedManager
    from .query_set import QuerySet
    from .record import Record
    from .ulabel import ULabel


# does not inherit from TracksRun because the Transform
# is needed to define a run
class Transform(SQLRecord, IsVersioned):
    """Data transformations such as scripts, notebooks, functions, or pipelines.

    A `transform` can be a function, a script, a notebook, or a
    pipeline. If you execute a transform, you generate a run
    (:class:`~lamindb.Run`). A run has inputs and outputs.

    Pipelines are typically created with a workflow manager (Nextflow, Snakemake,
    Prefect, Flyte, Dagster, redun, Airflow, ...).

    Transforms are versioned so that a given transform version maps on a given
    source code version.

    .. dropdown:: Can I sync transforms to git?

        If you set the environment variable `LAMINDB_SYNC_GIT_REPO` or set
        `ln.settings.sync_git_repo`, a script-like transform is
        synced to its hashed state in a git repository upon calling `ln.track()`::

            ln.settings.sync_git_repo = "https://github.com/laminlabs/lamindb"
            ln.track()

        If the hash isn't found in the git repository, an error is thrown.

        You can also create transforms that map pipelines via `Transform.from_git()`.

    The definition of transforms and runs is consistent with the OpenLineage
    specification where a `transform` would be called a "job" and a `run` a "run".

    Args:
        key: `str | None = None` A short name or path-like semantic key.
        kind: `TransformKind | None = "pipeline"` See :class:`~lamindb.base.types.TransformKind`.
        version: `str | None = None` A version string.
        description: `str | None = None` A description.
        reference: `str | None = None` A reference, e.g., a URL.
        reference_type: `str | None = None` A reference type, e.g., 'url'.
        source_code: `str | None = None` Source code of the transform.
        revises: `Transform | None = None` An old version of the transform.
        skip_hash_lookup: `bool = False` Skip the hash lookup so that a new transform is created even if a transform with the same hash already exists.

    See Also:
        :func:`~lamindb.track`
            Track a script or notebook run.
        :class:`~lamindb.Run`
            Executions of transforms.

    Notes:
        - :doc:`docs:track`
        - :doc:`docs:redun`
        - :doc:`docs:nextflow`
        - :doc:`docs:snakemake`

    Examples:

        Create a transform by running `ln.track()` in a notebook or a script::

            ln.track()

        Create a transform for a standalone function that acts as its own workflow::

            @ln.flow()
            def my_workflow():
                print("Hello, world!")

        Create a transform for a step in a workflow::

            @ln.step()
            def my_step():
                print("One step!")

        Create a transform for a pipeline::

            transform = ln.Transform(key="Cell Ranger", version="7.2.0", kind="pipeline").save()

        Create a transform by saving a Python or shell script or a notebook via the CLI::

            lamin save my_script.py
            lamin save my_script.sh
            lamin save my_notebook.ipynb

    """

    class Meta(SQLRecord.Meta, IsVersioned.Meta):
        abstract = False
        app_label = "lamindb"
        unique_together = ("key", "hash")

    _len_stem_uid: int = 12
    _len_full_uid: int = 16
    _name_field: str = "key"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=_len_full_uid
    )
    """Universal id."""
    # the max length equals the max length of an S3 key & the artifact key
    key: str = CharField(db_index=True, max_length=1024)
    """A name or "/"-separated path-like string.

    All transforms with the same key are part of the same version family.
    """
    # db_index on description because sometimes we query for equality in the case of artifacts
    description: str | None = TextField(null=True, db_index=True)
    """A description."""
    kind: TransformKind = CharField(
        max_length=20,
        db_index=True,
        default="pipeline",
    )
    """A string indicating the kind of transform (default `"pipeline"`).

    One of `"pipeline"`, `"notebook"`, `"script"`, or `"function"`.
    """
    source_code: str | None = TextField(null=True)
    """Source code of the transform."""
    hash: str | None = CharField(max_length=HASH_LENGTH, db_index=True, null=True)
    """Hash of the source code."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """Reference for the transform, e.g., a URL."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Reference type of the transform, e.g., 'url'."""
    environment: Artifact | None = models.ForeignKey(
        "Artifact", CASCADE, null=True, related_name="_environment_of_transforms"
    )
    """An environment for executing the transform."""
    plan: Artifact | None = models.ForeignKey(
        "Artifact",
        CASCADE,
        null=True,
        related_name="_plan_for_transforms",
        default=None,
    )
    """Optional plan artifact (e.g. default agent plan for this transform)."""
    runs: RelatedManager[Run]
    """Runs of this transform ← :attr:`~lamindb.Run.transform`."""
    ulabels: RelatedManager[ULabel] = models.ManyToManyField(
        "ULabel", through="TransformULabel", related_name="transforms"
    )
    """ULabel annotations of this transform ← :attr:`~lamindb.ULabel.transforms`."""
    linked_in_records: RelatedManager[Record] = models.ManyToManyField(
        "Record", through="RecordTransform", related_name="linked_transforms"
    )
    """This transform is linked in these records as a value ← :attr:`~lamindb.Record.linked_transforms`."""
    records: RelatedManager[Record]
    """Records that annotate this transform ← :attr:`~lamindb.Record.transforms`."""
    predecessors: RelatedManager[Transform] = models.ManyToManyField(
        "self",
        through="TransformTransform",
        symmetrical=False,
        related_name="successors",
    )
    """Preceding transforms ← :attr:`~lamindb.Transform.successors`."""
    successors: RelatedManager[Transform]
    """Subsequent transforms ← :attr:`~lamindb.Transform.predecessors`.

    Allows defining succeeding transforms. Is *not* necessary for data lineage, which is tracked automatically
    whenever an artifact or collection serves as an input for a run.
    """
    projects: RelatedManager[Project]
    """Linked projects ← :attr:`~lamindb.Project.transforms`."""
    references: RelatedManager[Reference]
    """Linked references ← :attr:`~lamindb.Reference.transforms`."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    updated_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of last update to record."""
    created_by: User = ForeignKey(
        User, PROTECT, default=current_user_id, related_name="created_transforms"
    )
    """Creator of record ← :attr:`~lamindb.User.created_transforms`."""
    ablocks: RelatedManager[TransformBlock]
    """Attached blocks ← :attr:`~lamindb.TransformBlock.transform`."""

    @overload
    def __init__(
        self,
        key: str | None = None,
        kind: TransformKind | None = None,
        version: str | None = None,
        description: str | None = None,
        reference: str | None = None,
        reference_type: str | None = None,
        source_code: str | None = None,
        revises: Transform | None = None,
        skip_hash_lookup: bool = False,
    ): ...

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
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a Transform"
            )
        key: str | None = kwargs.pop("key", None)
        description: str | None = kwargs.pop("description", None)
        revises: Transform | None = kwargs.pop("revises", None)
        version_tag: str | None = kwargs.pop("version_tag", kwargs.pop("version", None))
        kind: TransformKind | None = kwargs.pop("kind", None)
        type: TransformKind | None = kwargs.pop("type", None)
        if type is not None:
            warnings.warn(
                "`type` argument of transform was renamed to `kind` and will be removed in a future release.",
                DeprecationWarning,
                stacklevel=2,
            )
        kind = kind if kind is not None else (type if type is not None else "pipeline")
        reference: str | None = kwargs.pop("reference", None)
        reference_type: str | None = kwargs.pop("reference_type", None)
        branch = kwargs.pop("branch", None)
        branch_id = kwargs.pop("branch_id", 1)
        space = kwargs.pop("space", None)
        space_id = kwargs.pop("space_id", 1)
        skip_hash_lookup: bool = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        # below is internal use that we'll hopefully be able to eliminate
        uid: str | None = kwargs.pop("uid") if "uid" in kwargs else None
        source_code: str | None = (
            kwargs.pop("source_code") if "source_code" in kwargs else None
        )
        if not len(kwargs) == 0:
            raise ValueError(
                "Only key, description, version, kind, type, revises, reference, "
                f"reference_type can be passed, but you passed: {kwargs}"
            )
        if revises is None:
            # need to check uid before checking key
            if uid is not None:
                revises = (
                    Transform.objects.using(using_key)
                    .filter(uid__startswith=uid[:-4], is_latest=True)
                    .order_by("-created_at")
                    .first()
                )
            elif key is not None:
                candidate_for_revises = (
                    Transform.objects.using(using_key)
                    .filter(~Q(branch_id=-1), key=key, is_latest=True)
                    .order_by("-created_at")
                    .first()
                )
                if candidate_for_revises is not None:
                    revises = candidate_for_revises
                    if candidate_for_revises.source_code is None:
                        # no source code was yet saved, return the same transform
                        logger.important(
                            "no source code was yet saved, returning existing transform with same key"
                        )
                        uid = revises.uid
        if revises is not None and uid is not None and uid == revises.uid:
            if revises.key != key:
                logger.warning("ignoring inconsistent key")
            init_self_from_db(self, revises)
            update_attributes(self, {"description": description})
            return None
        if revises is not None and key is not None and revises.key != key:
            logger.important(f"renaming transform {revises.key} to {key}")
        new_uid, version_tag, key, description, revises = process_revises(
            revises, version_tag, key, description, Transform
        )
        # this is only because the user-facing constructor allows passing a uid
        # most others don't
        if uid is None:
            has_consciously_provided_uid = False
            uid = new_uid
        else:
            has_consciously_provided_uid = True
        hash = None
        if source_code is not None and not skip_hash_lookup:
            hash = hash_string(source_code)
            transform_candidate = Transform.objects.filter(
                ~Q(branch_id=-1),
                hash=hash,
                is_latest=True,
            ).first()
            if transform_candidate is not None:
                init_self_from_db(self, transform_candidate)
                update_attributes(self, {"description": description})
                if key is not None and transform_candidate.key != key:
                    logger.warning(
                        f"key {self.key} on existing transform differs from passed key {key}, keeping original key; update manually if needed or pass skip_hash_lookup if you want to duplicate the transform"
                    )
                return None
        super().__init__(  # type: ignore
            uid=uid,
            description=description,
            key=key,
            kind=kind,
            version_tag=version_tag,
            reference=reference,
            reference_type=reference_type,
            source_code=source_code,
            hash=hash,
            _has_consciously_provided_uid=has_consciously_provided_uid,
            revises=revises,
            branch=branch,
            branch_id=branch_id,
            space=space,
            space_id=space_id,
        )

    @classmethod
    def from_git(
        cls,
        url: str,
        path: str,
        key: str | None = None,
        version: str | None = None,
        entrypoint: str | None = None,
        branch: str | None = None,
        description: str | None = None,
        skip_hash_lookup: bool = False,
    ) -> Transform:
        """Create a transform from a path in a git repository.

        Args:
            url: URL of the git repository.
            path: Path to the file within the repository.
            key: Optional key for the transform.
            version: Optional version tag to checkout in the repository.
            entrypoint: One or several optional comma-separated entrypoints for the transform.
            branch: Optional branch to checkout.
            description: Optional description for the transform.
            skip_hash_lookup: Skip the hash lookup so that a new transform is created even if a transform with the same hash already exists.

        Examples:

            Create from a Nextflow repo and auto-infer the commit hash from its latest version::

                transform = ln.Transform.from_git(
                    url="https://github.com/openproblems-bio/task_batch_integration",
                    path="main.nf"
                ).save()

            Create from a Nextflow repo and checkout a specific version::

                transform = ln.Transform.from_git(
                    url="https://github.com/openproblems-bio/task_batch_integration",
                    path="main.nf",
                    version="v2.0.0"
                ).save()
                assert transform.version_tag == "v2.0.0"

            Create a *sliding transform* from a Nextflow repo's `dev` branch.
            Unlike a regular transform, a sliding transform doesn't pin a specific source code state,
            but adapts to whatever the referenced state on the branch is::

                transform = ln.Transform.from_git(
                    url="https://github.com/openproblems-bio/task_batch_integration",
                    path="main.nf",
                    branch="dev",
                    version="dev",
                ).save()

        Notes:

            A regular transform pins a specific source code state through its commit hash::

                transform.source_code
                #> repo: https://github.com/openproblems-bio/task_batch_integration
                #> path: main.nf
                #> commit: 68eb2ecc52990617dbb6d1bb5c7158d9893796bb

            A sliding transform infers the source code state from a branch::

                transform.source_code
                #> repo: https://github.com/openproblems-bio/task_batch_integration
                #> path: main.nf
                #> branch: dev

            If an entrypoint is provided, it is added to the source code below the path, e.g.::

                transform.source_code
                #> repo: https://github.com/openproblems-bio/task_batch_integration
                #> path: main.nf
                #> entrypoint: myentrypoint
                #> commit: 68eb2ecc52990617dbb6d1bb5c7158d9893796bb

            Note that you can pass a comma-separated list of entrypoints to the `entrypoint` argument.

        """
        from ..core._sync_git import get_and_validate_git_metadata

        url, commit_hash = get_and_validate_git_metadata(url, path, version, branch)
        if key is None:
            key = (
                url.split("/")[-2]
                + "/"
                + url.split("/")[-1].replace(".git", "")
                + "/"
                + path
            )
            logger.important(f"inferred key '{key}' from url & path")
        source_code = f"repo: {url}\npath: {path}"
        if entrypoint is not None:
            source_code += f"\nentrypoint: {entrypoint}"
        if branch is not None and version == branch:
            from urllib.parse import quote

            # sliding transform, no defined source code state
            source_code += f"\nbranch: {branch}"
            reference, reference_type = (
                f"{url}/tree/{quote(branch, safe='')}/{path}",
                "url",
            )
        else:
            # regular transform, defined source code state
            source_code += f"\ncommit: {commit_hash}"
            reference, reference_type = f"{url}/blob/{commit_hash}/{path}", "url"
        return Transform(
            key=key,
            kind="pipeline",
            version=version,
            description=description,
            reference=reference,
            reference_type=reference_type,
            source_code=source_code,
            skip_hash_lookup=skip_hash_lookup,
        )

    @property
    def latest_run(self) -> Run:
        """The latest run of this transform."""
        return self.runs.order_by("-started_at").first()

    @property
    @deprecated(new_name="kind")
    def type(self) -> TransformKind:
        return self.kind

    @type.setter
    def type(self, value: TransformKind):
        self.kind = value

    def view_lineage(self, with_successors: bool = False, distance: int = 5):
        """View lineage of transforms.

        Note that this only accounts for manually defined predecessors and successors.

        Auto-generate lineage through inputs and outputs of runs is not included.
        """
        from .has_parents import view_parents

        return view_parents(
            record=self,
            field="key",
            with_children=with_successors,
            distance=distance,
            attr_name="predecessors",
        )

    def _update_source_code_from_path(self, source_code_path: Path) -> None | str:
        _, transform_hash, _ = hash_file(source_code_path)  # ignore hash_type for now
        if self.hash is not None:
            # check if the hash of the transform source code matches
            if transform_hash != self.hash:
                response = input(
                    f"You are about to overwrite existing source code (hash '{self.hash}') for Transform('{self.uid}')."
                    f" Proceed? (y/n) "
                )
                if response == "y":
                    self.source_code = source_code_path.read_text()
                    self.hash = transform_hash
                else:
                    logger.warning("Please re-run `ln.track()` to make a new version")
                    return "rerun-the-notebook"
            else:
                logger.debug("source code is already saved")
        else:
            self.source_code = source_code_path.read_text()
            self.hash = transform_hash
        return None


def _permanent_delete_transforms(transforms: Transform | QuerySet) -> None:
    """Execute bulk DELETE on transforms (runs, then transforms). Used by QuerySet and single-transform paths."""
    from django.db.models import QuerySet as DjangoQuerySet

    from .project import TransformProject

    if isinstance(transforms, Transform):
        db = transforms._state.db or "default"
        qs = Transform.objects.using(db).filter(pk=transforms.pk)
    else:
        db = transforms.db or "default"
        qs = transforms
    objects = list(qs)
    if not objects:
        return
    _adjust_is_latest_when_deleting_is_versioned(objects)
    transform_ids = [o.pk for o in objects]
    TransformProject.objects.using(db).filter(transform_id__in=transform_ids).delete()
    Run.objects.using(db).filter(transform_id__in=transform_ids).delete(permanent=True)
    DjangoQuerySet.delete(qs)


class TransformTransform(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    successor: Transform = ForeignKey(
        "Transform", CASCADE, related_name="links_predecessor"
    )
    predecessor: Transform = ForeignKey(
        "Transform", CASCADE, related_name="links_successor"
    )
    config: dict | None = models.JSONField(default=None, null=True)
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now()
    )
    created_by: User = ForeignKey(
        "lamindb.User", PROTECT, default=current_user_id, related_name="+"
    )

    class Meta:
        app_label = "lamindb"
        unique_together = ("successor", "predecessor")
