from __future__ import annotations

from typing import TYPE_CHECKING, overload

import numpy as np
from django.db import models
from django.db.models import (
    CASCADE,
    PROTECT,
)
from lamindb_setup import _check_instance_setup

from lamindb.base import deprecated
from lamindb.base.fields import (
    BooleanField,
    CharField,
    DateTimeField,
    ForeignKey,
)
from lamindb.base.users import current_user_id
from lamindb.errors import InvalidArgument

from ..base.ids import base62_16
from .can_curate import CanCurate
from .sqlrecord import BaseSQLRecord, IsLink, SQLRecord

if TYPE_CHECKING:
    from datetime import datetime

    from .artifact import Artifact
    from .collection import Collection
    from .feature import FeatureValue
    from .project import Project
    from .query_set import QuerySet
    from .transform import Transform
    from .ulabel import ULabel


_TRACKING_READY: bool | None = None


class FeatureManager:
    """Feature manager."""

    pass


class FeatureManagerRun(FeatureManager):
    """Feature manager."""

    pass


def current_run() -> Run | None:
    global _TRACKING_READY

    if not _TRACKING_READY:
        _TRACKING_READY = _check_instance_setup()
    if _TRACKING_READY:
        import lamindb

        # also see get_run() in core._data
        run = lamindb._tracked.get_current_tracked_run()
        if run is None:
            run = lamindb.context.run
        return run
    else:
        return None


class TracksRun(models.Model):
    """Base class tracking latest run, creating user, and `created_at` timestamp."""

    class Meta:
        abstract = True

    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "lamindb.User",
        PROTECT,
        editable=False,
        default=current_user_id,
        related_name="+",
    )
    """Creator of record."""
    run: Run | None = ForeignKey(
        "lamindb.Run", PROTECT, null=True, default=current_run, related_name="+"
    )
    """Run that created record."""

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
        super().__init__(*args, **kwargs)


class TracksUpdates(models.Model):
    """Base class tracking previous runs and `updated_at` timestamp."""

    class Meta:
        abstract = True

    updated_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of last update to record."""

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
        super().__init__(*args, **kwargs)


class User(BaseSQLRecord, CanCurate):
    """Users.

    Every :class:`~lamindb.models.SQLRecord` has a `created_by` field that links to the creating user.

    All data in this registry is synchronized from LaminHub to ensure a universal
    user identity.

    Examples:

        Query a user by handle:

        >>> user = ln.User.get(handle="testuser1")
        >>> user
    """

    _name_field: str = "handle"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(editable=False, unique=True, db_index=True, max_length=8)
    """Universal id, valid across DB instances."""
    handle: str = CharField(max_length=30, unique=True, db_index=True)
    """User handle, valid across DB instances (required)."""
    name: str | None = CharField(max_length=150, db_index=True, null=True)
    """Full name (optional)."""  # has to match hub specification, where it's also optional
    created_artifacts: Artifact
    """Artifacts created by user."""
    created_transforms: Transform
    """Transforms created by user."""
    created_runs: Run
    """Runs created by user."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    updated_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of last update to record."""

    @overload
    def __init__(
        self,
        handle: str,
        email: str,
        name: str | None,
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
        super().__init__(*args, **kwargs)


class Run(SQLRecord):
    """Runs of transforms such as the execution of a script.

    A registry to store runs of transforms, such as an executation of a script.

    Args:
        transform: `Transform` A :class:`~lamindb.Transform` record.
        reference: `str | None = None` For instance, an external ID or a download URL.
        reference_type: `str | None = None` For instance, `redun_id`, `nextflow_id` or `url`.

    See Also:
        :meth:`~lamindb.core.Context.track`
            Track global runs & transforms for a notebook or script.

    Examples:

        Create a run record:

        >>> ln.Transform(key="Cell Ranger", version="7.2.0", type="pipeline").save()
        >>> transform = ln.Transform.get(key="Cell Ranger", version="7.2.0")
        >>> run = ln.Run(transform)

        Create a global run context for a custom transform:

        >>> ln.track(transform=transform)
        >>> ln.context.run  # globally available run

        Track a global run context for a notebook or script:

        >>> ln.track()  # Jupyter notebook metadata is automatically parsed
        >>> ln.context.run
    """

    _name_field: str = "started_at"

    features: FeatureManager = FeatureManagerRun  # type: ignore
    """Features manager.

    Run parameters are tracked via the `Feature` registry, just like all other variables.

    Guide: :ref:`track-run-parameters`

    Example::

        run.features.add_values({
            "learning_rate": 0.01,
            "input_dir": "s3://my-bucket/mydataset",
            "downsample": True,
            "preprocess_params": {
                "normalization_type": "cool",
                "subset_highlyvariable": True,
            },
        })
    """

    id: int = models.BigAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    # default uid was changed from base62_20 to base62_16 in 1.6.0
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=20, default=base62_16
    )
    """Universal id, valid across DB instances."""
    name: str | None = CharField(max_length=150, null=True)
    """A name."""
    transform: Transform = ForeignKey("Transform", CASCADE, related_name="runs")
    """The transform :class:`~lamindb.Transform` that is being run."""
    started_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Start time of run."""
    finished_at: datetime | None = DateTimeField(db_index=True, null=True, default=None)
    """Finished time of run."""
    # we don't want to make below a OneToOne because there could be the same trivial report
    # generated for many different runs
    report: Artifact | None = ForeignKey(
        "Artifact", PROTECT, null=True, related_name="_report_of", default=None
    )
    """Report of run, e.g.. n html file."""
    _logfile: Artifact | None = ForeignKey(
        "Artifact", PROTECT, null=True, related_name="_logfile_of", default=None
    )
    """Report of run, e.g.. n html file."""
    environment: Artifact | None = ForeignKey(
        "Artifact", PROTECT, null=True, related_name="_environment_of", default=None
    )
    """Computational environment for the run.

    For instance, `Dockerfile`, `docker image`, `requirements.txt`, `environment.yml`, etc.
    """
    input_artifacts: Artifact
    """The artifacts serving as input for this run.

    Related accessor: :attr:`~lamindb.Artifact.input_of_runs`.
    """
    output_artifacts: Artifact
    """The artifacts generated by this run.

    Related accessor: via :attr:`~lamindb.Artifact.run`
    """
    input_collections: Collection
    """The collections serving as input for this run."""
    output_collections: Collection
    """The collections generated by this run."""
    """Parameter values."""
    _feature_values: FeatureValue = models.ManyToManyField(
        "FeatureValue", through="RunFeatureValue", related_name="runs"
    )
    """Feature values."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """A reference like a URL or external ID (such as from a workflow manager)."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Type of reference such as a workflow manager execution ID."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of first creation. Mismatches ``started_at`` if the run is re-run."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=current_user_id, related_name="created_runs"
    )
    """Creator of run."""
    ulabels: ULabel = models.ManyToManyField(
        "ULabel", through="RunULabel", related_name="runs"
    )
    """ULabel annotations of this transform."""
    initiated_by_run: Run | None = ForeignKey(
        "Run", CASCADE, null=True, related_name="initiated_runs", default=None
    )
    """The run that triggered the current run.

    This is not a preceding run. The preceding runs ("predecessors") is the set
    of runs that produced the output artifacts that serve as the inputs for the
    present run.
    """
    initiated_runs: Run
    """Runs that were initiated by this run."""
    projects: Project
    """Linked projects."""
    _is_consecutive: bool | None = BooleanField(null=True)
    """Indicates whether code was consecutively executed. Is relevant for notebooks."""
    _status_code: int = models.SmallIntegerField(default=0, db_index=True)
    """Status code of the run.

    - 0: scheduled
    - 1: started
    - 2: errored
    - 3: aborted
    - 4: completed
    """

    @overload
    def __init__(
        self,
        transform: Transform,
        reference: str | None = None,
        reference_type: str | None = None,
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
        self.features = FeatureManager(self)  # type: ignore
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        # now we proceed with the user-facing constructor
        if len(args) > 1:
            raise ValueError("Only one non-keyword arg allowed: transform")
        transform: Transform = None
        if "transform" in kwargs or len(args) == 1:
            transform = kwargs.pop("transform") if len(args) == 0 else args[0]
        reference: str | None = kwargs.pop("reference", None)
        reference_type: str | None = kwargs.pop("reference_type", None)
        initiated_by_run: Run | None = kwargs.pop("initiated_by_run", None)
        if transform is None:
            raise TypeError("Pass transform parameter")
        if transform._state.adding:
            raise ValueError("Please save transform record before creating a run")

        super().__init__(  # type: ignore
            transform=transform,
            reference=reference,
            initiated_by_run=initiated_by_run,
            reference_type=reference_type,
        )

    def delete(self) -> None:
        """Delete."""
        delete_run_artifacts(self)
        super().delete()

    @property
    @deprecated("features")
    def params(self) -> FeatureManager:
        return self.features

    @classmethod
    def filter(
        cls,
        *queries,
        **expressions,
    ) -> QuerySet:
        """Query a set of artifacts.

        Args:
            *queries: `Q` expressions.
            **expressions: Params, fields, and values passed via the Django query syntax.

        See Also:
            - Guide: :doc:`docs:registries`

        Examples:

            Query by fields::

                ln.Run.filter(key="examples/my_file.parquet")

            Query by params::

                ln.Run.filter(hyperparam_x=100)
        """
        from ._feature_manager import filter_base
        from .feature import Feature
        from .query_set import QuerySet

        if expressions:
            keys_normalized = [key.split("__")[0] for key in expressions]
            field_or_feature_or_param = keys_normalized[0].split("__")[0]
            if field_or_feature_or_param in Run.__get_available_fields__():
                return QuerySet(model=cls).filter(*queries, **expressions)
            elif all(
                params_validated := Feature.validate(
                    keys_normalized, field="name", mute=True
                )
            ):
                return filter_base(FeatureManagerRun, **expressions)
            else:
                params = ", ".join(sorted(np.array(keys_normalized)[~params_validated]))
                message = f"feature names: {params}"
                fields = ", ".join(sorted(cls.__get_available_fields__()))
                raise InvalidArgument(
                    f"You can query either by available fields: {fields}\n"
                    f"Or fix invalid {message}"
                )
        else:
            return QuerySet(model=cls).filter(*queries, **expressions)


def delete_run_artifacts(run: Run) -> None:
    environment = None
    if run.environment is not None:
        environment = run.environment
        run.environment = None
    report = None
    if run.report is not None:
        report = run.report
        run.report = None
    if environment is not None or report is not None:
        run.save()
    if environment is not None:
        # only delete if there are no other runs attached to this environment
        if environment._environment_of.count() == 0:
            environment.delete(permanent=True)
    if report is not None:
        # only delete if there are no other runs attached to this environment
        if report._report_of.count() == 0:
            report.delete(permanent=True)


class RunFeatureValue(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    run: Run = ForeignKey(Run, CASCADE, related_name="links_featurevalue")
    # we follow the lower() case convention rather than snake case for link models
    featurevalue: FeatureValue = ForeignKey(
        "FeatureValue", PROTECT, related_name="links_run"
    )
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "lamindb.User", PROTECT, default=current_user_id, related_name="+"
    )
    """Creator of record."""

    class Meta:
        unique_together = ("run", "featurevalue")
