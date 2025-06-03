from __future__ import annotations

from typing import TYPE_CHECKING, Any, overload

from django.db import models
from django.db.models import CASCADE, PROTECT

from lamindb.base.fields import (
    BooleanField,
    CharField,
    ForeignKey,
    JSONField,
)
from lamindb.errors import FieldValidationError

from ..base.ids import base62_12, base62_16
from .artifact import Artifact
from .can_curate import CanCurate
from .feature import Feature
from .run import Run, TracksRun, TracksUpdates
from .sqlrecord import BaseSQLRecord, IsLink, SQLRecord, _get_record_kwargs
from .ulabel import ULabel

if TYPE_CHECKING:
    from .project import Project
    from .schema import Schema


class Record(SQLRecord, CanCurate, TracksRun, TracksUpdates):
    """Flexible records to register, e.g., samples, donors, cells, compounds, sequences.

    This is currently more convenient to use through the UI.

    A `Record` has a flexible schema: it can store data for arbitrary features.
    Changing the fields of a :class:`~lamindb.models.SQLRecord`, you need to modify the columns of the underlying table in the database.

    Args:
        name: `str` A name.
        description: `str` A description.

    See Also:
        :meth:`~lamindb.Sheet`
            Sheets to group records.
        :meth:`~lamindb.Feature`
            Dimensions of measurement.
        :attr:`~lamindb.Artifact.features`
            Feature manager for an artifact.
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
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type.

    For example, if a record "Compound" is a `type`, the actual compounds "darerinib", "tramerinib", would be instances of that `type`.
    """
    # naming convention in analogy with Schema
    components: Record = models.ManyToManyField(
        "Record", through="RecordRecord", symmetrical=False, related_name="composites"
    )
    """Record-like components of this record."""
    composites: Record
    """Record-like composites of this record."""
    sheet: Sheet | None = ForeignKey(
        "Sheet", CASCADE, null=True, related_name="records"
    )
    """Group records by sheet."""
    description: str | None = CharField(null=True)
    """A description (optional)."""
    artifacts: Artifact = models.ManyToManyField(
        Artifact, through="RecordArtifact", related_name="records"
    )
    """Linked artifacts."""
    runs: Run = models.ManyToManyField(Run, through="RecordRun", related_name="records")
    """Linked runs."""
    ulabels: ULabel = models.ManyToManyField(
        ULabel,
        through="RecordULabel",
        related_name="_records",  # in transition period
    )
    """Linked runs."""
    projects: Project
    """Linked projects."""

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
        sheet: Sheet = kwargs.pop("sheet", None)
        description: str | None = kwargs.pop("description", None)
        _skip_validation = kwargs.pop(
            "_skip_validation", True
        )  # should not validate records
        _aux = kwargs.pop("_aux", None)
        if len(kwargs) > 0:
            valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Record)])
            raise FieldValidationError(
                f"Only {valid_keywords} are valid keyword arguments"
            )
        super().__init__(
            name=name,
            type=type,
            is_type=is_type,
            sheet=sheet,
            description=description,
            _skip_validation=_skip_validation,
            _aux=_aux,
        )


class Sheet(SQLRecord, TracksRun, TracksUpdates):
    """Sheets to group records."""

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    id: int = models.AutoField(primary_key=True)
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=12, default=base62_12
    )
    """A universal random id, valid across DB instances."""
    name: str = CharField(db_index=True)
    """Name or title of sheet."""
    schema: Schema | None = ForeignKey(
        "Schema", CASCADE, null=True, related_name="sheets"
    )
    """A schema to enforce for the sheet (optional)."""
    description: str | None = CharField(null=True, db_index=True)
    """A description (optional)."""
    projects: Project
    """Linked projects."""

    @overload
    def __init__(
        self,
        name: str,
        schema: Schema | None = None,
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
        schema: Schema | None = kwargs.pop("schema", None)
        description: str | None = kwargs.pop("description", None)
        _skip_validation = kwargs.pop("_skip_validation", True)
        _aux = kwargs.pop("_aux", None)
        if len(kwargs) > 0:
            valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Record)])
            raise FieldValidationError(
                f"Only {valid_keywords} are valid keyword arguments"
            )
        super().__init__(
            name=name,
            schema=schema,
            description=description,
            _skip_validation=_skip_validation,
            _aux=_aux,
        )


class RecordJson(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_json")
    feature: Feature = ForeignKey(Feature, CASCADE, related_name="links_recordjson")
    value: Any = JSONField(default=None, db_default=None)

    class Meta:
        unique_together = ("record", "feature")


class RecordRecord(SQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(
        Record, CASCADE, related_name="values_record"
    )  # composite
    feature: Feature = ForeignKey(Feature, CASCADE, related_name="links_recordrecord")
    value: Record = ForeignKey(
        Record, PROTECT, related_name="links_record"
    )  # component

    class Meta:
        unique_together = ("record", "feature")


class RecordULabel(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_ulabel")
    feature: Feature = ForeignKey(Feature, CASCADE, related_name="links_recordulabel")
    value: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_record")

    class Meta:
        # allows linking exactly one record to one ulabel per feature, because we likely don't want to have Many
        unique_together = ("record", "feature")


class RecordRun(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_run")
    feature: Feature = ForeignKey(Feature, CASCADE, related_name="links_recordrun")
    value: Run = ForeignKey(Run, PROTECT, related_name="links_record")

    class Meta:
        # allows linking several records to a single run for the same feature because we'll likely need this
        unique_together = ("record", "feature")


class RecordArtifact(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_artifact")
    feature: Feature = ForeignKey(Feature, CASCADE, related_name="links_recordartifact")
    value: Artifact = ForeignKey(Artifact, PROTECT, related_name="links_record")

    class Meta:
        # allows linking several records to a single artifact for the same feature because we'll likely need this
        unique_together = ("record", "feature", "value")
