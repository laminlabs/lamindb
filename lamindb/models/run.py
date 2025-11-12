from __future__ import annotations

from typing import TYPE_CHECKING, overload

from django.db import models
from django.db.models import (
    CASCADE,
    PROTECT,
)
from lamin_utils import logger
from lamindb_setup import _check_instance_setup

from lamindb.base.fields import (
    BooleanField,
    CharField,
    DateTimeField,
    ForeignKey,
)
from lamindb.base.users import current_user_id

from ..base.ids import base62_16
from .can_curate import CanCurate
from .sqlrecord import BaseSQLRecord, IsLink, SQLRecord

if TYPE_CHECKING:
    from datetime import datetime

    from ._feature_manager import FeatureManager
    from .artifact import Artifact
    from .block import RunBlock
    from .collection import Collection
    from .feature import FeatureValue
    from .project import Project
    from .query_set import QuerySet
    from .record import Record
    from .transform import Transform
    from .ulabel import ULabel


_TRACKING_READY: bool | None = None


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

    This registry is automatically populated with user identities from LaminHub in case the user authenticates.

    Examples:

        Query a user by handle::

            user = ln.User.get(handle="testuser1")
    """

    class Meta:
        app_label = "lamindb"

    _name_field: str = "handle"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(editable=False, unique=True, db_index=True, max_length=8)
    """Universal id, valid across DB instances."""
    handle: str = CharField(max_length=30, unique=True, db_index=True)
    """User handle, valid across DB instances (required)."""
    name: str | None = CharField(max_length=150, db_index=True, null=True)
    """Full name (optional)."""  # has to match hub specification, where it's also optional
    linked_in_records: Record = models.ManyToManyField(
        "Record", through="RecordUser", related_name="linked_users"
    )
    """Records linked in this user."""
    artifacts: Artifact = models.ManyToManyField(
        "Artifact",
        through="ArtifactUser",
        through_fields=("user", "artifact"),
        related_name="users",
    )
    """Artifacts annotated with this user."""
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
        name: `str | None = None` An optional name.
        params: `dict | None = None` A dictionary of parameters.
        reference: `str | None = None` For instance, an external ID or a download URL.
        reference_type: `str | None = None` For instance, `redun_id`, `nextflow_id` or `url`.
        initiated_by_run: `Run | None = None` The run that triggers this run.

    See Also:
        :func:`~lamindb.track`
            Globally track a script or notebook run.
        :func:`~lamindb.tracked`
            Track a function with this decorator.

    Examples:

        Create a run record::

            ln.Transform(key="Cell Ranger", version="7.2.0", type="pipeline").save()
            transform = ln.Transform.get(key="Cell Ranger", version="7.2.0")
            run = ln.Run(transform)

        Create a global run context for a custom transform::

            ln.track(transform=transform)
            ln.context.run  # global run

        Track a global run context for a notebook or script::

            ln.track()
            ln.context.run  # global run

        You can pass parameters to `Run(transform, params=params)` or add them later::

            run.params = {
                "learning_rate": 0.01,
                "input_dir": "s3://my-bucket/mydataset",
                "downsample": True,
                "preprocess_params": {
                    "normalization_type": "cool",
                    "subset_highlyvariable": True,
                },
            }
            run.save()

        In contrast to `.params`, features are indexed in the `Feature` registry and can reference relational categorical values.
        If you want to link feature values, use::

            run.features.add_values({
                "experiment": "My experiment 1",
            })

        Guide: :ref:`track-run-parameters`
    """

    class Meta:
        app_label = "lamindb"

    _name_field: str = "started_at"
    _aux_fields: dict[str, tuple[str, type]] = {
        "0": ("cli_args", str),
    }

    id: int = models.BigAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    # default uid was changed from base62_20 to base62_16 in 1.6.0
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=20, default=base62_16
    )
    """Universal id, valid across DB instances."""
    name: str | None = CharField(max_length=150, null=True, db_index=True)
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
    input_records: Record
    """The collections serving as input for this run."""
    output_records: Record
    """The collections generated by this run."""
    params: dict = models.JSONField(null=True)
    """JSON-like parameters."""
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
    blocks: RunBlock
    """Blocks that annotate this run."""
    records: Record
    """Records that annotate this run."""
    linked_in_records: Record = models.ManyToManyField(
        "Record", through="RecordRun", related_name="linked_runs"
    )
    """This run is linked in these records as a value."""
    _is_consecutive: bool | None = BooleanField(null=True)
    """Indicates whether code was consecutively executed. Is relevant for notebooks."""
    _status_code: int = models.SmallIntegerField(
        default=-3, db_default=-3, db_index=True, null=True
    )
    """Status code of the run. See the status property for mapping to string."""

    @overload
    def __init__(
        self,
        transform: Transform,
        name: str | None = None,
        params: dict | None = None,
        reference: str | None = None,
        reference_type: str | None = None,
        initiated_by_run: Run | None = None,
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
        # now we proceed with the user-facing constructor
        if len(args) > 1:
            raise ValueError("Only one non-keyword arg allowed: transform")
        transform: Transform = None
        if "transform" in kwargs or len(args) == 1:
            transform = kwargs.pop("transform") if len(args) == 0 else args[0]
        name: str | None = kwargs.pop("name", None)
        params: dict | None = kwargs.pop("params", None)
        reference: str | None = kwargs.pop("reference", None)
        reference_type: str | None = kwargs.pop("reference_type", None)
        initiated_by_run: Run | None = kwargs.pop("initiated_by_run", None)
        if transform is None:
            raise TypeError("Pass transform parameter")
        if transform._state.adding:
            raise ValueError("Please save transform record before creating a run")
        if not len(kwargs) == 0:
            raise ValueError(
                f"Only transform, name, params, reference, reference_type, initiated_by_run can be passed, but you passed: {kwargs}"
            )
        super().__init__(  # type: ignore
            transform=transform,
            name=name,
            params=params,
            reference=reference,
            reference_type=reference_type,
            initiated_by_run=initiated_by_run,
        )

    @property
    def status(self) -> str:
        """Get status of run.

        Returns the status as a string, one of: `scheduled`, `re-started`, `started`, `completed`, `errored`, or `aborted`.

        Examples:

            See the status of a run::

                run.status
                #> 'completed'

        """
        if self._status_code is None:
            return "unknown"
        status_dict = {
            -3: "scheduled",
            -2: "re-started",
            -1: "started",
            0: "completed",
            1: "errored",
            2: "aborted",
        }
        return status_dict.get(self._status_code, "unknown")

    @property
    def cli_args(self) -> str | None:
        """CLI arguments if the run was invoked from the command line."""
        if self._aux is not None and "af" in self._aux and "0" in self._aux["af"]:  # type: ignore
            return self._aux["af"]["0"]  # type: ignore
        else:
            return None

    @cli_args.setter
    def cli_args(self, value: str) -> None:
        if not isinstance(value, str) or not (value):
            logger.warning("did not set empty or non-string cli_args")
            return
        self._aux = self._aux or {}
        self._aux.setdefault("af", {})["0"] = value

    @property
    def features(self) -> FeatureManager:
        """Manage annotations with features."""
        from ._feature_manager import FeatureManager

        return FeatureManager(self)

    def describe(self, return_str: bool = False) -> None | str:
        """Describe record including relations.

        Args:
            return_str: Return a string instead of printing.
        """
        from ._describe import describe_postgres_sqlite

        return describe_postgres_sqlite(self, return_str=return_str)

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
        # from Registry metaclass
        return type(cls).filter(cls, *queries, **expressions)


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
        # only delete if there are no other runs attached to this report
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
        app_label = "lamindb"
        unique_together = ("run", "featurevalue")
