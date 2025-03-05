from __future__ import annotations

from typing import TYPE_CHECKING, Any, overload

from django.db import models
from django.db.models import (
    CASCADE,
    PROTECT,
    Q,
)
from django.db.utils import IntegrityError
from lamindb_setup import _check_instance_setup
from lamindb_setup.core.hashing import HASH_LENGTH, hash_dict

from lamindb.base.fields import (
    BooleanField,
    CharField,
    DateTimeField,
    ForeignKey,
)
from lamindb.base.users import current_user_id
from lamindb.errors import ValidationError

from ..base.ids import base62_20
from .can_curate import CanCurate
from .record import BasicRecord, LinkORM, Record

if TYPE_CHECKING:
    from datetime import datetime

    from .artifact import Artifact
    from .collection import Collection
    from .project import Project
    from .schema import Schema
    from .transform import Transform
    from .ulabel import ULabel


_TRACKING_READY: bool | None = None


class ParamManager:
    """Param manager."""

    pass


class ParamManagerRun(ParamManager):
    """Param manager."""

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


class User(BasicRecord, CanCurate):
    """Users.

    All data in this registry is synced from `lamin.ai` to ensure a universal
    user identity. There is no need to manually create records.

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


class Param(Record, CanCurate, TracksRun, TracksUpdates):
    """Parameters of runs & models."""

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "name"

    name: str = CharField(max_length=100, db_index=True)
    dtype: str | None = CharField(db_index=True, null=True)
    """Data type ("num", "cat", "int", "float", "bool", "datetime").

    For categorical types, can define from which registry values are
    sampled, e.g., `cat[ULabel]` or `cat[bionty.CellType]`.
    """
    type: Param | None = ForeignKey("self", PROTECT, null=True, related_name="records")
    """Type of param (e.g., 'Pipeline', 'ModelTraining', 'PostProcessing').

    Allows to group features by type, e.g., all read outs, all metrics, etc.
    """
    records: Param
    """Records of this type."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    _expect_many: bool = models.BooleanField(default=False, db_default=False)
    """Indicates whether values for this param are expected to occur a single or multiple times for an artifact/run (default `False`).

    - if it's `False` (default), the values mean artifact/run-level values and a dtype of `datetime` means `datetime`
    - if it's `True`, the values are from an aggregation, which this seems like an edge case but when characterizing a model ensemble trained with different parameters it could be relevant
    """
    schemas: Schema = models.ManyToManyField(
        "Schema", through="SchemaParam", related_name="params"
    )
    """Feature sets linked to this feature."""
    # backward fields
    values: ParamValue
    """Values for this parameter."""

    def __init__(self, *args, **kwargs):
        from .feature import process_init_feature_param

        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None

        dtype = kwargs.get("dtype", None)
        kwargs = process_init_feature_param(args, kwargs, is_param=True)
        super().__init__(*args, **kwargs)
        dtype_str = kwargs.pop("dtype", None)
        if not self._state.adding:
            if not (
                self.dtype.startswith("cat")
                if dtype == "cat"
                else self.dtype == dtype_str
            ):
                raise ValidationError(
                    f"Feature {self.name} already exists with dtype {self.dtype}, you passed {dtype_str}"
                )


# FeatureValue behaves in many ways like a link in a LinkORM
# in particular, we don't want a _public field on it
# Also, we don't inherit from TracksRun because a ParamValue
# is typically created before a run is created and we want to
# avoid delete cycles (for Model params though it might be helpful)
class ParamValue(Record):
    """Parameter values.

    Is largely analogous to `FeatureValue`.
    """

    # we do not have a unique constraint on param & value because it leads to hashing errors
    # for large dictionaries: https://lamin.ai/laminlabs/lamindata/transform/jgTrkoeuxAfs0000
    # we do not hash values because we have `get_or_create` logic all over the place
    # and also for checking whether the (param, value) combination exists
    # there does not seem an issue with querying for a dict-like value
    # https://lamin.ai/laminlabs/lamindata/transform/jgTrkoeuxAfs0001
    _name_field: str = "value"

    param: Param = ForeignKey(Param, CASCADE, related_name="values")
    """The dimension metadata."""
    value: Any = (
        models.JSONField()
    )  # stores float, integer, boolean, datetime or dictionaries
    """The JSON-like value."""
    # it'd be confusing and hard to populate a run here because these
    # values are typically created upon creating a run
    # hence, ParamValue does _not_ inherit from TracksRun but manually
    # adds created_at & created_by
    # because ParamValue cannot be updated, we don't need updated_at
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        User, PROTECT, default=current_user_id, related_name="+"
    )
    """Creator of record."""
    hash: str = CharField(max_length=HASH_LENGTH, null=True, db_index=True)

    class Meta:
        constraints = [
            # For simple types, use direct value comparison
            models.UniqueConstraint(
                fields=["param", "value"],
                name="unique_simple_param_value",
                condition=Q(hash__isnull=True),
            ),
            # For complex types (dictionaries), use hash
            models.UniqueConstraint(
                fields=["param", "hash"],
                name="unique_complex_param_value",
                condition=Q(hash__isnull=False),
            ),
        ]

    @classmethod
    def get_or_create(cls, param, value):
        # Simple types: int, float, str, bool
        if isinstance(value, (int, float, str, bool)):
            try:
                return cls.objects.create(param=param, value=value, hash=None), False
            except IntegrityError:
                return cls.objects.get(param=param, value=value), True

        # Complex types: dict, list
        else:
            hash = hash_dict(value)
            try:
                return cls.objects.create(param=param, value=value, hash=hash), False
            except IntegrityError:
                return cls.objects.get(param=param, hash=hash), True


class Run(Record):
    """Runs.

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

    params: ParamManager = ParamManagerRun  # type: ignore
    """Param manager.

    Guide: :ref:`track-run-parameters`

    Example::

        run.params.add_values({
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
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=20, default=base62_20
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
    _param_values: ParamValue = models.ManyToManyField(
        ParamValue, through="RunParamValue", related_name="runs"
    )
    """Parameter values."""
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
        self.params = ParamManager(self)  # type: ignore
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
        report.delete(permanent=True)


class RunParamValue(BasicRecord, LinkORM):
    id: int = models.BigAutoField(primary_key=True)
    run: Run = ForeignKey(Run, CASCADE, related_name="+")
    # we follow the lower() case convention rather than snake case for link models
    paramvalue: ParamValue = ForeignKey(ParamValue, PROTECT, related_name="+")
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "lamindb.User", PROTECT, default=current_user_id, related_name="+"
    )
    """Creator of record."""

    class Meta:
        unique_together = ("run", "paramvalue")
