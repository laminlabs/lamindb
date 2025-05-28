from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, overload

from django.db import models
from django.db.models import PROTECT
from lamin_utils import logger
from lamindb_setup.core.hashing import HASH_LENGTH, hash_string

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

    from .project import Project, Reference
    from .ulabel import ULabel


# does not inherit from TracksRun because the Transform
# is needed to define a run
class Transform(SQLRecord, IsVersioned):
    """Data transformations such as scripts, notebooks, functions, or pipelines.

    A "transform" can refer to a Python function, a script, a notebook, or a
    pipeline. If you execute a transform, you generate a run
    (:class:`~lamindb.Run`). A run has inputs and outputs.

    A pipeline is typically created with a workflow tool (Nextflow, Snakemake,
    Prefect, Flyte, MetaFlow, redun, Airflow, ...) and stored in a versioned
    repository.

    Transforms are versioned so that a given transform version maps on a given
    source code version.

    .. dropdown:: Can I sync transforms to git?

        If you switch on
        :attr:`~lamindb.core.Settings.sync_git_repo` a script-like transform is
        synched to its hashed state in a git repository upon calling `ln.track()`.

        >>> ln.settings.sync_git_repo = "https://github.com/laminlabs/lamindb"
        >>> ln.track()

    The definition of transforms and runs is consistent the OpenLineage
    specification where a :class:`~lamindb.Transform` record would be called a
    "job" and a :class:`~lamindb.Run` record a "run".

    Args:
        name: `str` A name or title.
        key: `str | None = None` A short name or path-like semantic key.
        type: `TransformType | None = "pipeline"` See :class:`~lamindb.base.types.TransformType`.
        revises: `Transform | None = None` An old version of the transform.

    See Also:
        :meth:`~lamindb.core.Context.track`
            Globally track a script, notebook or pipeline run.
        :class:`~lamindb.Run`
            Executions of transforms.

    Notes:
        - :doc:`docs:track`
        - :doc:`docs:data-flow`
        - :doc:`docs:redun`
        - :doc:`docs:nextflow`
        - :doc:`docs:snakemake`

    Examples:

        Create a transform for a pipeline:

        >>> transform = ln.Transform(key="Cell Ranger", version="7.2.0", type="pipeline").save()

        Create a transform from a notebook:

        >>> ln.track()

        View predecessors of a transform:

        >>> transform.view_lineage()
    """

    class Meta(SQLRecord.Meta, IsVersioned.Meta):
        abstract = False
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
    key: str | None = CharField(db_index=True, null=True)
    """A name or "/"-separated path-like string.

    All transforms with the same key are part of the same version family.
    """
    description: str | None = CharField(db_index=True, null=True)
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

    @overload
    def __init__(
        self,
        name: str,
        key: str | None = None,
        type: TransformType | None = None,
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
        key: str | None = kwargs.pop("key", None)
        description: str | None = kwargs.pop("description", None)
        revises: Transform | None = kwargs.pop("revises", None)
        version: str | None = kwargs.pop("version", None)
        type: TransformType | None = kwargs.pop("type", "pipeline")
        reference: str | None = kwargs.pop("reference", None)
        reference_type: str | None = kwargs.pop("reference_type", None)
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
                    .filter(key=key, is_latest=True)
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
            transform_candidate = Transform.filter(
                hash=hash, is_latest=True
            ).one_or_none()
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
        )

    @property
    def name(self) -> str:
        """Name of the transform.

        Splits `key` on `/` and returns the last element.
        """
        return self.key.split("/")[-1]

    @property
    def latest_run(self) -> Run:
        """The latest run of this transform."""
        return self.runs.order_by("-started_at").first()

    def delete(self) -> None:
        """Delete."""
        # query all runs and delete their artifacts
        runs = Run.filter(transform=self)
        for run in runs:
            delete_run_artifacts(run)
        # at this point, all artifacts have been taken care of
        # we can now leverage CASCADE delete
        super().delete()

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
