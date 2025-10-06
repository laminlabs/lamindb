from __future__ import annotations

from typing import TYPE_CHECKING, Any, overload

from django.db import models
from django.db.models import CASCADE, PROTECT
from lamin_utils import logger
from lamindb_setup.core import deprecated

from lamindb.base.fields import (
    BooleanField,
    CharField,
    ForeignKey,
    JSONField,
    TextField,
)
from lamindb.errors import FieldValidationError

from ..base.ids import base62_16
from .artifact import Artifact
from .can_curate import CanCurate
from .feature import Feature
from .has_parents import HasParents, _query_relatives
from .query_set import reorder_subset_columns_in_df
from .run import Run, TracksRun, TracksUpdates, User, current_run
from .sqlrecord import BaseSQLRecord, IsLink, SQLRecord, _get_record_kwargs
from .transform import Transform
from .ulabel import ULabel

if TYPE_CHECKING:
    import pandas as pd

    from .blocks import RunBlock
    from .project import Project, Reference
    from .query_set import QuerySet
    from .schema import Schema


class Record(SQLRecord, CanCurate, TracksRun, TracksUpdates, HasParents):
    """Metadata records for labeling and organizing entities in sheets.

    Is useful to manage samples, donors, cells, compounds, sequences.

    Args:
        name: `str` A name.
        description: `str` A description.
        type: `Record | None = None` The type of this record.
        is_type: `bool = False` Whether this record is a type (a record that
            classifies other records).
        schema: `Schema | None = None` A schema to enforce for a type (optional).
        reference: `str | None = None` For instance, an external ID or a URL.
        reference_type: `str | None = None` For instance, `"url"`.

    See Also:
        :meth:`~lamindb.Feature`
            Dimensions of measurement (e.g. column of a sheet, attribute of a record).

    Examples:

        Create a record type and then instances of that type::

            sample_type = Record(name="Sample", is_type=True).save()
            sample1 = Record(name="Sample 1", type=sample_type).save()
            sample2 = Record(name="Sample 2", type=sample_type).save()

        You can then annotate artifacts and other entities with these records, e.g.::

            artifact.records.add(sample1)

        To query artifacts by records::

            ln.Artifact.filter(records=sample1).to_dataframe()

        Through the UI can assign attributes to records in form of features. The Python API also allows to
        assign features programmatically, but is currently still low-level::

            feature = ln.Feature(name="age", type="int").save()
            sample1.values_record.create(feature=feature, value=42)
            sample2.values_record.create(feature=feature, value=23)

        Records can also model flexible ontologies through their parents-children relationships::

            cell_type = Record(name="CellType", is_type=True).save()
            t_cell = Record(name="T Cell", type=cell_type).save()
            cd4_t_cell = Record(name="CD4+ T Cell", type=cell_type).save()
            t_cell.children.add(cd4_t_cell)

        Often, a label is measured *within* a dataset. For instance, an artifact
        might characterize 2 species of the Iris flower (`"setosa"` &
        `"versicolor"`) measured by a `"species"` feature. For such cases, you can use
        :class:`~lamindb.curators.DataFrameCurator` to automatically parse, validate, and
        annotate with labels that are contained in `DataFrame` objects.

    .. note::

        If you work with complex entities like cell lines, cell types, tissues,
        etc., consider using the pre-defined biological registries in
        :mod:`bionty` to label artifacts & collections.

        If you work with biological samples, likely, the only sustainable way of
        tracking metadata, is to create a custom schema module.

    .. note::

        A `Record` has a flexible schema: it can store data for arbitrary features.
        By contrast, if you want to change the fields of a :class:`~lamindb.models.SQLRecord`, you need to modify the columns of the underlying table in the database.
        The latter is more efficient for large datasets and you can customize it through modules like the `bionty` or `wetlab` module.

    """

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False
        app_label = "lamindb"
        constraints = [
            models.UniqueConstraint(
                fields=["name", "type", "space"], name="unique_name_type_space"
            )
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
    """Records of this type (can only be non-empty if `is_type` is `True`)."""
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
    """A schema to enforce for a type (optional).

    This is mostly parallel to the `schema` attribute of `Artifact`.

    If `is_type` is `True`, the schema is used to enforce certain features for each records of this type.
    """
    # naming convention in analogy to Schema
    components: Record = models.ManyToManyField(
        "Record", through="RecordRecord", symmetrical=False, related_name="composites"
    )
    """Record-like components of this record."""
    composites: Record
    """Record-like composites of this record."""
    parents: ULabel = models.ManyToManyField(
        "self", symmetrical=False, related_name="children"
    )
    """Parent entities of this record.

    For advanced use cases, you can build an ontology under a given `type`.

    Say, if you modeled `CellType` as a `Record`, you would introduce a type `CellType` and model the hiearchy of cell types under it.
    """
    children: ULabel
    """Child entities of this record.

    Reverse accessor for parents.
    """
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
    projects: Project
    """Projects that annotate this record."""
    references: Reference
    """References that annotate this record."""
    values_json: RecordJson
    """JSON values (for lists, dicts, etc.)."""
    values_record: RecordRecord
    """Record values with their features."""
    values_ulabel: RecordULabel
    """ULabel values with their features."""
    values_user: RecordUser
    """User values with their features."""
    values_run: RecordRun
    """Run values with their features."""
    values_artifact: RecordArtifact
    """Artifact values with their features."""
    values_reference: Reference
    """Reference values with their features."""
    values_project: Project
    """Project values with their features."""
    linked_runs: Run = models.ManyToManyField(
        Run, through="RecordRun", related_name="records"
    )
    """Runs linked in this record as values."""
    linked_users: User = models.ManyToManyField(
        User, through="RecordUser", related_name="records"
    )
    """Users linked in this record as values."""
    linked_ulabels: ULabel = models.ManyToManyField(
        ULabel,
        through="RecordULabel",
        related_name="linked_in_records",
    )
    """ULabels linked in this record as values."""
    linked_artifacts: Artifact = models.ManyToManyField(
        Artifact, through="RecordArtifact", related_name="linked_in_records"
    )
    """Artifacts linked in this record as values."""
    linked_projects: Project
    """Projects linked in this record as values."""
    linked_references: Reference
    """References linked in this record as values."""
    blocks: RunBlock
    """Blocks that annotate this record."""

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
    def is_form(self) -> bool:
        """Check if record is a form (a record type with a validating schema)."""
        return self.schema is not None and self.is_type

    def query_children(self) -> QuerySet:
        """Query all children of a record.

        While `.children` retrieves the direct children, this method
        retrieves all descendants of a record type.
        """
        return _query_relatives([self], "children", self.__class__)  # type: ignore

    def query_records(self) -> QuerySet:
        """Query all records of a type.

        While `.records` retrieves the direct children, this method
        retrieves all descendants of a record type.
        """
        return _query_relatives([self], "records", self.__class__)  # type: ignore

    def type_to_dataframe(self) -> pd.DataFrame:
        """Export all instances of this record type to a pandas DataFrame."""
        assert self.is_type, "Only types can be exported as dataframes"  # noqa: S101
        df = self.query_records().to_dataframe(features="queryset")
        df.columns.values[0] = "__lamindb_record_uid__"
        df.columns.values[1] = "__lamindb_record_name__"
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
        format: dict[str, Any] = {"suffix": ".csv"} if key.endswith(".csv") else {}
        format["index"] = False
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
            format=format,
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


class RecordRun(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_run")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordrun")
    value: Run = ForeignKey(Run, PROTECT, related_name="links_record")

    class Meta:
        # allows linking several records to a single run for the same feature because we'll likely need this
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


class RecordArtifact(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_artifact")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordartifact")
    value: Artifact = ForeignKey(Artifact, PROTECT, related_name="links_in_record")

    class Meta:
        # allows linking several records to a single artifact for the same feature because we'll likely need this
        app_label = "lamindb"
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
        app_label = "lamindb"
        unique_together = ("artifact", "record", "feature")
