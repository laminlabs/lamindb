from __future__ import annotations

from typing import TYPE_CHECKING, Any, overload

from django.db import models
from django.db.models import CASCADE, PROTECT
from lamin_utils import logger

from lamindb.base.fields import (
    BooleanField,
    CharField,
    ForeignKey,
    JSONField,
)
from lamindb.errors import FieldValidationError

from ..base.ids import base62_16
from .artifact import Artifact
from .can_curate import CanCurate
from .feature import Feature
from .has_parents import _query_relatives
from .query_set import reorder_subset_columns_in_df
from .run import Run, TracksRun, TracksUpdates
from .sqlrecord import BaseSQLRecord, IsLink, SQLRecord, _get_record_kwargs
from .transform import Transform
from .ulabel import ULabel

if TYPE_CHECKING:
    import pandas as pd

    from .project import Person, Project, Reference
    from .query_set import QuerySet
    from .schema import Schema


class Record(SQLRecord, CanCurate, TracksRun, TracksUpdates):
    """Flexible records as you find them in Excel-like sheets.

    Useful register, e.g., samples, donors, cells, compounds, sequences.

    This is currently more convenient to use through the UI.

    A `Record` has a flexible schema: it can store data for arbitrary features.
    Changing the fields of a :class:`~lamindb.models.SQLRecord`, you need to modify the columns of the underlying table in the database.

    Args:
        name: `str` A name.
        description: `str` A description.

    See Also:
        :meth:`~lamindb.Feature`
            Dimensions of measurement (e.g. column of a sheet).
    """

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "name"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=16, default=base62_16
    )
    """A universal random id, valid across DB instances."""
    name: str = CharField(max_length=150, db_index=True, null=True)
    """Name or title of record (optional)."""
    type: Record | None = ForeignKey("self", PROTECT, null=True, related_name="records")
    """Type of record, e.g., `Sample`, `Donor`, `Cell`, `Compound`, `Sequence`.

    Allows to group records by type, e.g., all samples, all donors, all cells, all compounds, all sequences.
    """
    records: Record
    """Records of this type (can only be non-empty if `is_type` is `True`)."""
    is_type: bool = BooleanField(default=False, db_index=True)
    """Indicates if record is a `type`.

    For example, if a record "Compound" is a `type`, the actual compounds "darerinib", "tramerinib", would be instances of that `type`.
    """
    schema: Schema | None = ForeignKey(
        "Schema", CASCADE, null=True, related_name="records"
    )
    """A schema to enforce for a type (optional).

    This is mostly parallel to the `schema` attribute of `Artifact`.

    If `is_type` is `True`, the schema is used to enforce certain features for each records of this type.
    """
    # naming convention in analogy with Schema
    components: Record = models.ManyToManyField(
        "Record", through="RecordRecord", symmetrical=False, related_name="composites"
    )
    """Record-like components of this record."""
    composites: Record
    """Record-like composites of this record."""
    description: str | None = CharField(null=True)
    """A description (optional)."""
    linked_artifacts: Artifact = models.ManyToManyField(
        Artifact, through="RecordArtifact", related_name="linked_in_records"
    )
    """Linked artifacts."""
    artifacts: Artifact = models.ManyToManyField(
        Artifact, through="ArtifactRecord", related_name="records"
    )
    """Annotated artifacts."""
    linked_runs: Run = models.ManyToManyField(
        Run, through="RecordRun", related_name="records"
    )
    """Linked runs."""
    run: Run | None = ForeignKey(
        Run,
        PROTECT,
        related_name="output_records",
        null=True,
        default=None,
        editable=False,
    )
    """Run that created the record."""
    input_of_runs: Run = models.ManyToManyField(Run, related_name="input_records")
    """Runs that use this record as an input."""
    ulabels: ULabel = models.ManyToManyField(
        ULabel,
        through="RecordULabel",
        related_name="_records",  # in transition period
    )
    """Linked runs."""
    linked_projects: Project
    """Linked projects."""
    linked_references: Reference
    """Linked references."""
    linked_people: Person
    """Linked people."""

    @overload
    def __init__(
        self,
        name: str,
        type: Record | None = None,
        is_type: bool = False,
        description: str | None = None,
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
        schema = kwargs.pop("schema", None)
        branch = kwargs.pop("branch", None)
        branch_id = kwargs.pop("branch_id", 1)
        space = kwargs.pop("space", None)
        space_id = kwargs.pop("space_id", 1)
        _skip_validation = kwargs.pop(
            "_skip_validation", False
        )  # should not validate records
        _aux = kwargs.pop("_aux", None)
        if len(kwargs) > 0:
            valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Record)])
            raise FieldValidationError(
                f"Only {valid_keywords} are valid keyword arguments"
            )
        if schema and not is_type:
            logger.important("passing schema, treating as type")
            is_type = True
        super().__init__(
            name=name,
            type=type,
            is_type=is_type,
            description=description,
            schema=schema,
            branch=branch,
            branch_id=branch_id,
            space=space,
            space_id=space_id,
            _skip_validation=_skip_validation,
            _aux=_aux,
        )

    @property
    def is_form(self) -> bool:
        """Check if record is a form (a record type with a validating schema)."""
        return self.schema is not None and self.is_type

    def query_children(self) -> QuerySet:
        """Query all children of a record type recursively.

        While `.records` retrieves the direct children, this method
        retrieves all descendants of a record type.
        """
        return _query_relatives([self], "records", self.__class__)  # type: ignore

    def to_pandas(self) -> pd.DataFrame:
        """Export all children of a record type recursively to a pandas DataFrame."""
        assert self.is_type, "Only types can be exported as dataframes"  # noqa: S101
        df = self.query_children().df(features="queryset")
        df.columns.values[0] = "__lamindb_record_uid__"
        df.columns.values[1] = "__lamindb_record_name__"
        if self.schema is not None:
            desired_order = self.schema.members.list("name")  # only members is ordered!
        else:
            # sort alphabetically for now
            desired_order = df.columns[2:].tolist()
            desired_order.sort()
        df = reorder_subset_columns_in_df(df, desired_order, position=0)  # type: ignore
        return df.sort_index()  # order by id for now

    def to_artifact(self, key: str = None) -> Artifact:
        """Export all children of a record type as a `.csv` artifact."""
        from lamindb.core._context import context

        assert self.is_type, "Only types can be exported as artifacts"  # noqa: S101
        if key is None:
            file_suffix = ".csv"
            key = f"sheet_exports/{self.name}{file_suffix}"
        description = f": {self.description}" if self.description is not None else ""
        format: dict[str, Any] = {"suffix": ".csv"} if key.endswith(".csv") else {}
        format["index"] = False
        transform, _ = Transform.objects.get_or_create(
            key="__lamindb_record_export__", type="function"
        )
        run = Run(transform, initiated_by_run=context.run).save()
        run.input_records.add(self)
        return Artifact.from_df(
            self.to_pandas(),
            key=key,
            description=f"Export of sheet {self.uid}{description}",
            schema=self.schema,
            format=format,
            run=run,
        ).save()


class RecordJson(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_json")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordjson")
    value: Any = JSONField(default=None, db_default=None)

    class Meta:
        unique_together = ("record", "feature")  # a list is modeled as a list in json


class RecordRecord(SQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(
        Record, CASCADE, related_name="values_record"
    )  # composite
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordrecord")
    value: Record = ForeignKey(
        Record, PROTECT, related_name="links_record"
    )  # component

    class Meta:
        unique_together = ("record", "feature", "value")


class RecordULabel(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_ulabel")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordulabel")
    value: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_record")

    class Meta:
        # allows linking exactly one record to one ulabel per feature, because we likely don't want to have Many
        unique_together = ("record", "feature", "value")


class RecordRun(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_run")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordrun")
    value: Run = ForeignKey(Run, PROTECT, related_name="links_record")

    class Meta:
        # allows linking several records to a single run for the same feature because we'll likely need this
        unique_together = ("record", "feature", "value")


class RecordArtifact(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_artifact")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordartifact")
    value: Artifact = ForeignKey(Artifact, PROTECT, related_name="links_in_record")

    class Meta:
        # allows linking several records to a single artifact for the same feature because we'll likely need this
        unique_together = ("record", "feature", "value")


# like ArtifactULabel, for annotation
class ArtifactRecord(BaseSQLRecord, IsLink):
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
        unique_together = ("artifact", "record", "feature")
