from __future__ import annotations

from typing import TYPE_CHECKING, Any, overload

from django.db import models
from django.db.models import CASCADE, PROTECT
from lamin_utils import logger
from lamindb_setup.core import deprecated

from lamindb.base.fields import (
    BooleanField,
    CharField,
    DateTimeField,
    ForeignKey,
    JSONField,
    TextField,
)
from lamindb.errors import FieldValidationError

from ..base.ids import base62_16
from .artifact import Artifact
from .can_curate import CanCurate
from .feature import Feature
from .has_parents import HasParents, _query_ancestors_of_fk, _query_relatives
from .query_set import (
    QuerySet,
    SQLRecordList,
    encode_lamindb_fields_as_columns,
    get_default_branch_ids,
    reorder_subset_columns_in_df,
)
from .run import Run, TracksRun, TracksUpdates, User, current_run, current_user_id
from .sqlrecord import BaseSQLRecord, IsLink, SQLRecord, _get_record_kwargs
from .transform import Transform
from .ulabel import ULabel

if TYPE_CHECKING:
    from datetime import datetime

    import pandas as pd

    from ._feature_manager import FeatureManager
    from .block import RunBlock
    from .project import Project, RecordProject, RecordReference, Reference
    from .schema import Schema


class Record(SQLRecord, CanCurate, TracksRun, TracksUpdates, HasParents):
    """Flexible metadata records for labeling and organizing entities.

    Useful for managing samples, donors, cells, compounds, sequences, and other custom entities.

    Args:
        name: `str` A name.
        description: `str` A description.
        type: `Record | None = None` The type of this record.
        is_type: `bool = False` Whether this record is a type (a record that
            classifies other records).
        schema: `Schema | None = None` A schema defining allowed features for records of this type. Only applicable when `is_type=True`.
        reference: `str | None = None` For instance, an external ID or a URL.
        reference_type: `str | None = None` For instance, `"url"`.

    See Also:
        :meth:`~lamindb.Feature`
            Dimensions of measurement (e.g. column of a sheet, attribute of a record).

    Examples:

        Create a **record** and annotate an :class:`~lamindb.Artifact`::

            sample1 = ln.Record(name="Sample 1").save()
            artifact.records.add(sample1)

        Group several records under a **record type**::

            experiment_type = ln.Record(name="Experiment", is_type=True).save()
            experiment1 = ln.Record(name="Experiment 1", type=experiment_type).save()
            experiment2 = ln.Record(name="Experiment 2", type=experiment_type).save()

        Export all records of that type to dataframe::

            experiment_type.records.to_dataframe()
            #>              name   ...
            #>      Experiment 1   ...
            #>      Experiment 2   ...

        Add **features** to a record::

            gc_content = ln.Feature(name="gc_content", dtype=float).save()
            experiment = ln.Feature(name="experiment", dtype=experiment_type).save()
            sample1.features.add_values({
                "gc_content": 0.5,
                "experiment": "Experiment 1",
            })

        **Constrain features** by using a :class:`~lamindb.Schema`, creating a **sheet**::

            schema = ln.Schema([gc_content, experiment], name="sample_schema").save()
            sheet = ln.Record(name="Sample", is_type=True, schema=schema).save()  # add schema to type
            sample2 = ln.Record(name="Sample 2", type=sheet).save()
            sample2.features.add_values({"gc_content": 0.6})  # raises ValidationError because experiment is missing

        Query records by features::

            ln.Record.filter(gc_content=0.55)     # exact match
            ln.Record.filter(gc_content__gt=0.5)  # greater than
            ln.Record.filter(type=sheet)          # just the record on the sheet

        Model **custom ontologies** through their parents/children attributes::

            cell_type = ln.Record(name="CellType", is_type=True).save()
            t_cell = ln.Record(name="T Cell", type=cell_type).save()
            cd4_t_cell = ln.Record(name="CD4+ T Cell", type=cell_type).save()
            t_cell.children.add(cd4_t_cell)

        If you work with basic biological entities like cell lines, cell types, tissues,
        consider building on the public biological ontologies in :mod:`bionty`.

    .. dropdown:: What is the difference between `Record` and `SQLRecord`?

        The features of a `Record` are flexible: you can dynamically define features and add features to a record.
        The fields of a `SQLRecord` are fixed: you need to define them in code and then migrate the underlying database.

        You can configure a `SQLRecord` by subclassing it in a custom schema, for example, as done here: `github.com/laminlabs/wetlab <https://github.com/laminlabs/wetlab>`__

    """

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False
        app_label = "lamindb"
        constraints = [
            models.UniqueConstraint(
                fields=["name", "type", "space"],
                name="unique_record_name_type_space",
                condition=~models.Q(branch_id=-1),
            )
            # also see raw SQL constraints for `is_type` and `type` FK validity in migrations
        ]

    _name_field: str = "name"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=16, default=base62_16
    )
    """A universal random id, valid across DB instances."""
    name: str = CharField(max_length=150, db_index=True, null=True)
    """Name or title of record (optional).

    Names for a given `type` and `space` are constrained to be unique.
    """
    type: Record | None = ForeignKey("self", PROTECT, null=True, related_name="records")
    """Type of record, e.g., `Sample`, `Donor`, `Cell`, `Compound`, `Sequence`.

    Allows to group records by type, e.g., all samples, all donors, all cells, all compounds, all sequences.
    """
    records: Record
    """If a type (`is_type=True`), records of this `type`."""
    is_type: bool = BooleanField(default=False, db_index=True)
    """Indicates if record is a `type`.

    For example, if a record "Compound" is a `type`, the actual compounds "darerinib", "tramerinib", would be instances of that `type`.
    """
    description: str | None = TextField(null=True)
    """A description."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """A simple reference like a URL or external ID."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Type of simple reference."""
    schema: Schema | None = ForeignKey(
        "Schema", CASCADE, null=True, related_name="records"
    )
    """A schema to enforce for a type.

    This is analogous to the `schema` attribute of an `Artifact`.
    If `is_type` is `True`, the schema is used to enforce features for each record of this type.
    """
    # naming convention in analogy to Schema, but probably better would linked_records in analogy with other
    # record relationships
    components: Record = models.ManyToManyField(
        "Record", through="RecordRecord", symmetrical=False, related_name="composites"
    )
    """Records linked in this record as a value."""
    composites: Record  # consider renaming to linked_in_records in LaminDB 2
    """Records linking this record as a value. Is reverse accessor for `components`."""
    parents: Record = models.ManyToManyField(
        "self", symmetrical=False, related_name="children"
    )
    """Ontological parents of this record.

    You can build an ontology under a given `type`. For example, introduce a type `CellType` and model the hiearchy of cell types under it via `parents` and `children`.
    """
    children: Record
    """Ontological children of this record. Is reverse accessor for `parents`."""
    # this is handled manually here because we want to se the related_name attribute
    # (this doesn't happen via inheritance of TracksRun, everything else is the same)
    run: Run | None = ForeignKey(
        Run,
        PROTECT,
        related_name="output_records",
        null=True,
        default=current_run,
        editable=False,
    )
    """Run that created the record."""
    input_of_runs: Run = models.ManyToManyField(Run, related_name="input_records")
    """Runs that use this record as an input."""
    artifacts: Artifact = models.ManyToManyField(
        Artifact, through="ArtifactRecord", related_name="records"
    )
    """Artifacts annotated by this record."""
    runs: Run = models.ManyToManyField(Run, through="RunRecord", related_name="records")
    """Runs annotated by this record."""
    projects: Project
    """Projects that annotate this record."""
    references: Reference
    """References that annotate this record."""
    linked_runs: Run
    """Runs linked in this record as values."""
    linked_ulabels: ULabel
    """ULabels linked in this record as values."""
    linked_artifacts: Artifact
    """Artifacts linked in this record as values."""
    linked_projects: Project
    """Projects linked in this record as values."""
    linked_references: Reference
    """References linked in this record as values."""
    linked_users: User
    """Users linked in this record as values."""
    blocks: RunBlock
    """Blocks that annotate this record."""
    values_json: RecordJson
    """JSON values `(record_id, feature_id, value)`."""
    values_record: RecordRecord
    """Record values with their features `(record_id, feature_id, value_id)`."""
    values_ulabel: RecordULabel
    """ULabel values with their features `(record_id, feature_id, value_id)`."""
    values_user: RecordUser
    """User values with their features `(record_id, feature_id, value_id)`."""
    values_run: RecordRun
    """Run values with their features `(record_id, feature_id, value_id)`."""
    values_artifact: RecordArtifact
    """Artifact values with their features `(record_id, feature_id, value_id)`."""
    values_reference: RecordReference
    """Reference values with their features `(record_id, feature_id, value_id)`."""
    values_project: RecordProject
    """Project values with their features `(record_id, feature_id, value_id)`."""

    @overload
    def __init__(
        self,
        name: str,
        type: Record | None = None,
        is_type: bool = False,
        description: str | None = None,
        schema: Schema | None = None,
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
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if len(args) > 0:
            raise ValueError("Only one non-keyword arg allowed")
        name: str = kwargs.pop("name", None)
        type: str | None = kwargs.pop("type", None)
        is_type: bool = kwargs.pop("is_type", False)
        description: str | None = kwargs.pop("description", None)
        schema: Schema | None = kwargs.pop("schema", None)
        reference: str | None = kwargs.pop("reference", None)
        reference_type: str | None = kwargs.pop("reference_type", None)
        branch = kwargs.pop("branch", None)
        branch_id = kwargs.pop("branch_id", 1)
        space = kwargs.pop("space", None)
        space_id = kwargs.pop("space_id", 1)
        _skip_validation = kwargs.pop("_skip_validation", False)
        _aux = kwargs.pop("_aux", None)
        if len(kwargs) > 0:
            valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Record)])
            raise FieldValidationError(
                f"Only {valid_keywords} are valid keyword arguments"
            )
        if type and not type.is_type:
            raise ValueError(
                f"You can only assign a record of `is_type=True` as `type` to another record, but this doesn't have it: {type}"
            )
        if schema and not is_type:
            logger.important("passing schema, treating as type")
            is_type = True
        super().__init__(
            name=name,
            type=type,
            is_type=is_type,
            description=description,
            reference=reference,
            reference_type=reference_type,
            schema=schema,
            branch=branch,
            branch_id=branch_id,
            space=space,
            space_id=space_id,
            _skip_validation=_skip_validation,
            _aux=_aux,
        )

    @property
    def features(self) -> FeatureManager:
        """Manage annotations with features."""
        from ._feature_manager import FeatureManager

        return FeatureManager(self)

    @property
    def is_form(self) -> bool:
        """Check if record is a form, i.e., `self.is_type and self.schema is not None`."""
        return self.schema is not None and self.is_type

    def query_parents(self) -> QuerySet:
        """Query all parents of a record recursively.

        While `.parents` retrieves the direct parents, this method
        retrieves all ancestors of a record type.
        """
        return _query_relatives([self], "parents")  # type: ignore

    def query_children(self) -> QuerySet:
        """Query all children of a record recursively.

        While `.children` retrieves the direct children, this method
        retrieves all descendants of a record type.
        """
        return _query_relatives([self], "children")  # type: ignore

    def query_records(self) -> QuerySet:
        """Query all records of a type recursively.

        While `.records` retrieves the direct children, this method
        retrieves all descendants of a record type.
        """
        return _query_relatives([self], "records")  # type: ignore

    def query_types(self) -> SQLRecordList:
        """Query types of a record recursively.

        While `.type` retrieves the direct type, this method
        retrieves all ancestors of that `type`.
        """
        return _query_ancestors_of_fk(self, "type")  # type: ignore

    def type_to_dataframe(self, recurse: bool = False) -> pd.DataFrame:
        """Export all instances of this record type to a pandas DataFrame.

        This is almost equivalent to::

            ln.Record.filter(type=sample_type).to_dataframe(include="features")

        `type_to_dataframe()` ensures that the columns are ordered according to the schema of the type and encodes fields like `uid` and `name`.

        Args:
            recurse: `bool = False` Whether to include records of sub-types recursively.
        """
        assert self.is_type, "Only types can be exported as dataframes"  # noqa: S101

        branch_ids = get_default_branch_ids()
        qs = (
            self.query_records()
            if recurse
            else self.records.filter(branch_id__in=branch_ids)
        )
        df = qs.to_dataframe(features="queryset", order_by="id")
        encoded_id = encode_lamindb_fields_as_columns(self.__class__, "id")
        encoded_uid = encode_lamindb_fields_as_columns(self.__class__, "uid")
        encoded_name = encode_lamindb_fields_as_columns(self.__class__, "name")
        # encode the django id, uid and name fields
        if df.index.name == "id":
            df.index.name = encoded_id
        if "uid" in df.columns and encoded_uid not in df.columns:
            df = df.rename(columns={"uid": encoded_uid})
        if "name" in df.columns and encoded_name not in df.columns:
            df = df.rename(columns={"name": encoded_name})
        if self.schema is not None:
            desired_order = self.schema.members.to_list(
                "name"
            )  # only members is ordered!
        else:
            # sort alphabetically for now
            desired_order = df.columns[2:].tolist()
            desired_order.sort()
        df = reorder_subset_columns_in_df(df, desired_order, position=0)  # type: ignore
        return df.sort_index()  # order by id for now

    @deprecated("type_to_dataframe")
    def to_pandas(self) -> pd.DataFrame:
        return self.type_to_dataframe()

    def to_artifact(self, key: str = None) -> Artifact:
        """Calls `type_to_dataframe()` to create an artifact."""
        from lamindb.core._context import context

        assert self.is_type, "Only types can be exported as artifacts"  # noqa: S101
        if key is None:
            file_suffix = ".csv"
            key = f"sheet_exports/{self.name}{file_suffix}"
        description = f": {self.description}" if self.description is not None else ""
        transform, _ = Transform.objects.get_or_create(
            key="__lamindb_record_export__", type="function"
        )
        run = Run(transform, initiated_by_run=context.run).save()
        run.input_records.add(self)
        return Artifact.from_dataframe(
            self.type_to_dataframe(),
            key=key,
            description=f"Export of sheet {self.uid}{description}",
            schema=self.schema,
            csv_kwargs={"index": False},
            run=run,
        ).save()


class RecordJson(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_json")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordjson")
    value: Any = JSONField(default=None, db_default=None)

    class Meta:
        app_label = "lamindb"
        unique_together = ("record", "feature")  # a list is modeled as a list in json


class RecordRecord(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(
        Record, CASCADE, related_name="values_record"
    )  # composite
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordrecord")
    value: Record = ForeignKey(
        Record, PROTECT, related_name="links_record"
    )  # component

    class Meta:
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


class RecordULabel(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_ulabel")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordulabel")
    value: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_record")

    class Meta:
        # allows linking exactly one record to one ulabel per feature, because we likely don't want to have Many
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


class RecordUser(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_user")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recorduser")
    value: User = ForeignKey(User, PROTECT, related_name="links_record")

    class Meta:
        # allows linking exactly one record to one user per feature, because we likely don't want to have Many
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for storing run-like values in records
class RecordRun(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_run")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordrun")
    value: Run = ForeignKey(Run, PROTECT, related_name="links_in_record")

    class Meta:
        # allows linking several records to a single run for the same feature because we'll likely need this
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for annotating runs with records
class RunRecord(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    run: Run = ForeignKey(Run, CASCADE, related_name="links_record")
    record: Record = ForeignKey(Record, PROTECT, related_name="links_run")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_runrecord")
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    created_by: User = ForeignKey(
        "lamindb.User", PROTECT, default=current_user_id, related_name="+"
    )

    class Meta:
        app_label = "lamindb"
        unique_together = ("run", "record", "feature")


# for storing artifact-like values in records
class RecordArtifact(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_artifact")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordartifact")
    value: Artifact = ForeignKey(Artifact, PROTECT, related_name="links_in_record")

    class Meta:
        # allows linking several records to a single artifact for the same feature because we'll likely need this
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for annotating artifacts with records
class ArtifactRecord(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="links_record")
    record: Record = ForeignKey(Record, PROTECT, related_name="links_artifact")
    feature: Feature = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_artifactrecord"
    )
    label_ref_is_name: bool | None = BooleanField(null=True)
    feature_ref_is_name: bool | None = BooleanField(null=True)

    class Meta:
        # allows linking several records to a single artifact for the same feature because we'll likely need this
        app_label = "lamindb"
        unique_together = ("artifact", "record", "feature")
