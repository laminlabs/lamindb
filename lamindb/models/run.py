from __future__ import annotations

import os
import subprocess
import sys
from typing import TYPE_CHECKING, overload

from django.db import models
from django.db.models import (
    CASCADE,
    PROTECT,
)
from lamin_utils import logger
from lamindb_setup import _check_instance_setup
from lamindb_setup import settings as setup_settings

from lamindb.base.fields import (
    BooleanField,
    CharField,
    DateTimeField,
    ForeignKey,
)
from lamindb.base.users import current_user_id
from lamindb.base.utils import strict_classmethod

from ..base.uids import base62_16
from .can_curate import CanCurate
from .query_set import BasicQuerySet, QuerySet
from .sqlrecord import BaseSQLRecord, IsLink, SQLRecord

if TYPE_CHECKING:
    from datetime import datetime

    from ._feature_manager import FeatureManager
    from .artifact import Artifact
    from .block import RunBlock
    from .collection import Collection
    from .feature import JsonValue
    from .project import Project
    from .query_manager import RelatedManager
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
        run = lamindb.core._functions.get_current_tracked_run()
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
    linked_in_records: RelatedManager[Record] = models.ManyToManyField(
        "Record", through="RecordUser", related_name="linked_users"
    )
    """This user is linked in these records as a value."""
    artifacts: RelatedManager[Artifact] = models.ManyToManyField(
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
        uid: str,
        handle: str,
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


class Run(SQLRecord, TracksUpdates):
    """Runs of transforms such as the executions of a script.

    Args:
        transform: :class:`~lamindb.Transform` A data transformation object.
        name: `str | None = None` A name.
        params: `dict | None = None` A dictionary of parameters.
        reference: `str | None = None` For instance, an external ID or URL.
        reference_type: `str | None = None` For instance, `redun_id`, `nextflow_id` or `url`.
        initiated_by_run: `Run | None = None` The `run` that triggers this `run`.

    See Also:
        :func:`~lamindb.track`
            Globally track a script or notebook run.
        :func:`~lamindb.step`
            Track a function executionwith this decorator.

    Examples:

        Create a run record::

            ln.Transform(key="Cell Ranger", version="7.2.0", kind="pipeline").save()
            transform = ln.Transform.get(key="Cell Ranger", version="7.2.0")
            run = ln.Run(transform)

        Track a global run of a notebook or script::

            ln.track()
            ln.context.run  # global run object

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

    id: int = models.BigAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    # default uid was changed from base62_20 to base62_16 in 1.6.0
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=20, default=base62_16
    )
    """Universal id, valid across DB instances."""
    name: str | None = CharField(max_length=150, null=True, db_index=True)
    """An optional name for this run."""
    transform: Transform = ForeignKey("Transform", CASCADE, related_name="runs")
    """The transform that is being run ← :attr:`~lamindb.Transform.runs`."""
    entrypoint: str | None = CharField(max_length=255, null=True, db_index=True)
    """The entrypoint of the transform.

    This could be a function name or the entry point of a CLI or workflow manager.
    """
    started_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """The time this run started."""
    finished_at: datetime | None = DateTimeField(db_index=True, null=True, default=None)
    """The time this run finished or aborted."""
    # we don't want to make below a OneToOne because there could be the same trivial report
    # generated for many different runs
    report: Artifact | None = ForeignKey(
        "Artifact", PROTECT, null=True, related_name="_report_of", default=None
    )
    """The report of this run such as an `.html` or `.txt` file."""
    environment: Artifact | None = ForeignKey(
        "Artifact", PROTECT, null=True, related_name="_environment_of", default=None
    )
    """The computational environment for this run.

    For instance, `Dockerfile`, `docker image`, `requirements.txt`, `environment.yml`, etc.
    """
    input_artifacts: RelatedManager[Artifact]
    """The artifacts serving as input for this run, via :attr:`~lamindb.Artifact.input_of_runs`.
    """
    output_artifacts: Artifact
    """The artifacts generated by this run, via :attr:`~lamindb.Artifact.run`.
    """
    recreated_artifacts: RelatedManager[Artifact]
    """The output artifacts that were recreated by this run but originally created in another run, via :attr:`~lamindb.Artifact.recreating_runs`.

    Artifacts are considered recreated if they are reloaded due to a hash lookup match for an existing artifact.
    """
    input_records: RelatedManager[Record]
    """The collections serving as input for this run, via :attr:`~lamindb.Record.input_of_runs`."""
    output_records: Record
    """The collections generated by this run, via :attr:`~lamindb.Record.run`."""
    input_collections: RelatedManager[Collection]
    """The collections serving as input for this run, via :attr:`~lamindb.Collection.input_of_runs`."""
    output_collections: Collection
    """The collections generated by this run, via :attr:`~lamindb.Collection.run`."""
    recreated_collections: RelatedManager[Collection]
    """The output collections that were recreated by this run but originally created in another run, via :attr:`~lamindb.Collection.recreating_runs`.

    Artifacts are considered recreated if they are reloaded due to a hash lookup match for an existing artifact.
    """
    params: dict = models.JSONField(null=True)
    """Parameters (plain JSON values)."""
    json_values: RelatedManager[JsonValue] = models.ManyToManyField(
        "JsonValue", through="RunJsonValue", related_name="runs"
    )
    """Feature-indexed JSON values ← :attr:`~lamindb.JsonValue.runs`."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """A reference like a URL or an external ID such as from a workflow manager."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """The type of the `reference` such as a workflow manager execution ID."""
    cli_args: str | None = CharField(max_length=1024, null=True, default=None)
    """CLI arguments if the run was invoked from the command line."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """The time of creation of this run."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=current_user_id, related_name="created_runs"
    )
    """The creator of this run ← :attr:`~lamindb.User.created_runs`."""
    ulabels: RelatedManager[ULabel] = models.ManyToManyField(
        "ULabel", through="RunULabel", related_name="runs"
    )
    """The ulabels annotating this run ← :attr:`~lamindb.ULabel.runs`."""
    initiated_by_run: Run | None = ForeignKey(
        "Run", CASCADE, null=True, related_name="initiated_runs", default=None
    )
    """The run that initiated this run, via :attr:`~lamindb.Run.initiated_runs`."""
    initiated_runs: Run
    """The runs that were initiated by this run."""
    projects: RelatedManager[Project]
    """The projects annotating this run ← :attr:`~lamindb.Project.runs`."""
    ablocks: RunBlock
    """The blocks annotating this run ← :attr:`~lamindb.RunBlock.run`."""
    records: RelatedManager[Record]
    """The records annotating this run, via :attr:`~lamindb.Record.runs`."""
    linked_in_records: RelatedManager[Record] = models.ManyToManyField(
        "Record", through="RecordRun", related_name="linked_runs"
    )
    """This run is linked in these records as a value, via :attr:`~lamindb.Record.linked_runs`."""
    artifacts: RelatedManager[Artifact] = models.ManyToManyField(
        "Artifact", through="ArtifactRun", related_name="runs"
    )
    """The artifacts annotated by this run, via :attr:`~lamindb.Artifact.runs`."""
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
        entrypoint: str | None = None,
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
        entrypoint: str | None = kwargs.pop("entrypoint", None)
        params: dict | None = kwargs.pop("params", None)
        reference: str | None = kwargs.pop("reference", None)
        reference_type: str | None = kwargs.pop("reference_type", None)
        initiated_by_run: Run | None = kwargs.pop("initiated_by_run", None)
        report: Artifact | None = kwargs.pop("report", None)
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
            entrypoint=entrypoint,
            params=params,
            reference=reference,
            reference_type=reference_type,
            initiated_by_run=initiated_by_run,
            report=report,
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
    def features(self) -> FeatureManager:
        """Manage annotations with features."""
        from ._feature_manager import FeatureManager

        return FeatureManager(self)

    @strict_classmethod
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


def _permanent_delete_runs(runs: Run | QuerySet) -> None:
    """Execute bulk DELETE on runs and spawn artifact cleanup. Used by QuerySet and single-run paths."""
    if isinstance(runs, Run):
        db = runs._state.db or "default"
        first_run_uid = runs.uid
        artifact_ids = []
        if runs.environment_id:
            artifact_ids.append(runs.environment_id)
        if runs.report_id:
            artifact_ids.append(runs.report_id)
        super(BaseSQLRecord, runs).delete()
    else:
        db = runs.db or "default"
        rows = list(runs.values_list("uid", "report_id", "environment_id"))
        if rows:
            first_run_uid = rows[0][0]
        else:
            return
        artifact_ids = list({aid for r in rows for aid in r[1:3] if aid is not None})
        super(BasicQuerySet, runs).delete()
    if artifact_ids:
        ids_str = ",".join(map(str, artifact_ids))
        instance = db if db not in (None, "default") else setup_settings.instance.slug
        # spawn background subprocess to delete orphaned report/env artifacts
        cmd: list[str] = [
            sys.executable,
            "-m",
            "lamindb.models._run_cleanup",
            "--instance",
            instance,
            "--ids",
            ids_str,
            "--run-uid",
            first_run_uid,
        ]
        proc = subprocess.Popen(
            cmd,
            start_new_session=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=os.environ,
        )
        log_path = setup_settings.cache_dir / f"run_cleanup_logs_{first_run_uid}.txt"
        logger.debug(
            f"spawned run cleanup subprocess (pid={proc.pid}): {log_path}\n{' '.join(cmd)}"
        )


class RunJsonValue(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    run: Run = ForeignKey(Run, CASCADE, related_name="links_jsonvalue")
    # we follow the lower() case convention rather than snake case for link models
    jsonvalue: JsonValue = ForeignKey("JsonValue", PROTECT, related_name="links_run")
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
        unique_together = ("run", "jsonvalue")
