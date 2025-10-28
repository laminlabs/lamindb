from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, overload

from django.db import models
from django.db.models import PROTECT, Q
from lamin_utils import logger
from lamindb_setup.core.hashing import HASH_LENGTH, hash_string

from lamindb.base import deprecated
from lamindb.base.fields import (
    CharField,
    DateTimeField,
    ForeignKey,
    TextField,
)
from lamindb.base.users import current_user_id

from ..models._is_versioned import process_revises
from ._is_versioned import IsVersioned
from .run import Run, User, delete_run_artifacts
from .sqlrecord import SQLRecord, init_self_from_db, update_attributes

if TYPE_CHECKING:
    from datetime import datetime

    from lamindb.base.types import TransformType

    from .block import TransformBlock
    from .project import Project, Reference
    from .ulabel import ULabel


def delete_transform_relations(transform: Transform):
    from .project import TransformProject

    # query all runs and delete their associated report and env artifacts
    runs = Run.filter(transform=transform)
    for run in runs:
        delete_run_artifacts(run)
    # CASCADE doesn't do the job below because run_id might be protected through run__transform=self
    # hence, proactively delete the label links
    qs = TransformProject.filter(transform=transform)
    if qs.exists():
        qs.delete()
    # at this point, all artifacts have been taken care of
    # and one can now leverage CASCADE delete


# does not inherit from TracksRun because the Transform
# is needed to define a run
class Transform(SQLRecord, IsVersioned):
    """Data transformations such as scripts, notebooks, functions, or pipelines.

    A `transform` can be a function, a script, a notebook, or a
    pipeline. If you execute a transform, you generate a run
    (:class:`~lamindb.Run`). A run has inputs and outputs.

    Pipelines are typically created with a workflow tool (Nextflow, Snakemake,
    Prefect, Flyte, Dagster, redun, Airflow, ...).

    Transforms are versioned so that a given transform version maps on a given
    source code version.

    .. dropdown:: Can I sync transforms to git?

        If you switch on
        :attr:`~lamindb.core.Settings.sync_git_repo` a script-like transform is
        synched to its hashed state in a git repository upon calling `ln.track()`::

            ln.settings.sync_git_repo = "https://github.com/laminlabs/lamindb"
            ln.track()

        Alternatively, you create transforms that map pipelines via `Transform.from_git()`.

    The definition of transforms and runs is consistent with the OpenLineage
    specification where a `transform` would be called a "job" and a `run` a "run".

    Args:
        key: `str | None = None` A short name or path-like semantic key.
        type: `TransformType | None = "pipeline"` See :class:`~lamindb.base.types.TransformType`.
        version: `str | None = None` A version string.
        description: `str | None = None` A description.
        reference: `str | None = None` A reference, e.g., a URL.
        reference_type: `str | None = None` A reference type, e.g., 'url'.
        source_code: `str | None = None` Source code of the transform.
        revises: `Transform | None = None` An old version of the transform.

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

        Create a transform for a pipeline::

            transform = ln.Transform(key="Cell Ranger", version="7.2.0", type="pipeline").save()

        Create a transform from a notebook::

            ln.track()

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
    # the fact that key is nullable is consistent with Artifact
    # it might turn out that there will never really be a use case for this
    # but there likely also isn't much harm in it except for the mixed type
    # max length for key is 1014 and equals the max lenght of an S3 key & artifact key
    key: str | None = CharField(db_index=True, null=True, max_length=1024)
    """A name or "/"-separated path-like string.

    All transforms with the same key are part of the same version family.
    """
    # db_index on description because sometimes we query for equality in the case of artifacts
    description: str | None = TextField(null=True, db_index=True)
    """A description."""
    type: TransformType = CharField(
        max_length=20,
        db_index=True,
        default="pipeline",
    )
    """:class:`~lamindb.base.types.TransformType` (default `"pipeline"`)."""
    source_code: str | None = TextField(null=True)
    """Source code of the transform."""
    hash: str | None = CharField(max_length=HASH_LENGTH, db_index=True, null=True)
    """Hash of the source code."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """Reference for the transform, e.g., a URL."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Reference type of the transform, e.g., 'url'."""
    runs: Run
    """Runs of this transform."""
    ulabels: ULabel = models.ManyToManyField(
        "ULabel", through="TransformULabel", related_name="transforms"
    )
    """ULabel annotations of this transform."""
    predecessors: Transform = models.ManyToManyField(
        "self", symmetrical=False, related_name="successors"
    )
    """Preceding transforms.

    Allows to _manually_ define predecessors. Is typically not necessary as data lineage is
    automatically tracked via runs whenever an artifact or collection serves as an input for a run.
    """
    successors: Transform
    """Subsequent transforms.

    See :attr:`~lamindb.Transform.predecessors`.
    """
    projects: Project
    """Linked projects."""
    references: Reference
    """Linked references."""
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
    """Creator of record."""
    _template: Transform | None = ForeignKey(
        "Transform", PROTECT, related_name="_derived_from", default=None, null=True
    )
    """Creating template."""
    blocks: TransformBlock
    """Blocks that annotate this artifact."""

    @overload
    def __init__(
        self,
        key: str | None = None,
        type: TransformType | None = None,
        version: str | None = None,
        description: str | None = None,
        reference: str | None = None,
        reference_type: str | None = None,
        source_code: str | None = None,
        revises: Transform | None = None,
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
        version: str | None = kwargs.pop("version", None)
        type: TransformType | None = kwargs.pop("type", "pipeline")
        reference: str | None = kwargs.pop("reference", None)
        reference_type: str | None = kwargs.pop("reference_type", None)
        branch = kwargs.pop("branch", None)
        branch_id = kwargs.pop("branch_id", 1)
        space = kwargs.pop("space", None)
        space_id = kwargs.pop("space_id", 1)
        using_key = kwargs.pop("using_key", None)
        if "name" in kwargs:
            if key is None:
                key = kwargs.pop("name")
                warnings.warn(
                    f"`name` will be removed soon, please pass '{key}' to `key` instead",
                    FutureWarning,
                    stacklevel=2,
                )
            else:
                # description wasn't exist, so no check necessary
                description = kwargs.pop("name")
                warnings.warn(
                    f"`name` will be removed soon, please pass '{description}' to `description` instead",
                    FutureWarning,
                    stacklevel=2,
                )
        # below is internal use that we'll hopefully be able to eliminate
        uid: str | None = kwargs.pop("uid") if "uid" in kwargs else None
        source_code: str | None = (
            kwargs.pop("source_code") if "source_code" in kwargs else None
        )
        if not len(kwargs) == 0:
            raise ValueError(
                "Only key, description, version, type, revises, reference, "
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
        new_uid, version, key, description, revises = process_revises(
            revises, version, key, description, Transform
        )
        # this is only because the user-facing constructor allows passing a uid
        # most others don't
        if uid is None:
            has_consciously_provided_uid = False
            uid = new_uid
        else:
            has_consciously_provided_uid = True
        hash = None
        if source_code is not None:
            hash = hash_string(source_code)
            transform_candidate = Transform.objects.filter(
                ~Q(branch_id=-1),
                hash=hash,
                is_latest=True,
            ).first()
            if transform_candidate is not None:
                init_self_from_db(self, transform_candidate)
                update_attributes(self, {"key": key, "description": description})
                return None
        super().__init__(  # type: ignore
            uid=uid,
            description=description,
            key=key,
            type=type,
            version=version,
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

    @property
    @deprecated("key")
    def name(self) -> str:
        """Name of the transform.

        Splits `key` on `/` and returns the last element.
        """
        return self.key.split("/")[-1]

    def describe(self, return_str: bool = False) -> None | str:
        """Describe record including relations.

        Args:
            return_str: Return a string instead of printing.
        """
        from ._describe import describe_postgres_sqlite

        return describe_postgres_sqlite(self, return_str=return_str)

    @classmethod
    def from_git(
        cls,
        url: str,
        path: str,
        key: str | None = None,
        version: str | None = None,
        entrypoint: str | None = None,
        branch: str | None = None,
    ) -> Transform:
        """Create a transform from a path in a git repository.

        Args:
            url: URL of the git repository.
            path: Path to the file within the repository.
            key: Optional key for the transform.
            version: Optional version tag to checkout in the repository.
            entrypoint: Optional entrypoint for the transform.
            branch: Optional branch to checkout.

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
                assert transform.version == "v2.0.0"

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
            type="pipeline",
            version=version,
            reference=reference,
            reference_type=reference_type,
            source_code=source_code,
        )

    @property
    def latest_run(self) -> Run:
        """The latest run of this transform."""
        return self.runs.order_by("-started_at").first()

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
