from typing import TYPE_CHECKING, Optional

from django.db import models
from django.db.models import PROTECT, ManyToManyField

from lamindb.base.fields import (
    BooleanField,
    CharField,
    ForeignKey,
    IntegerField,
)

from .artifact import Artifact
from .feature import Feature
from .record import (
    BasicRecord,
    Record,
)
from .run import Param, TracksRun, TracksUpdates

if TYPE_CHECKING:
    from .project import Project

from django.core.exceptions import ValidationError

from ..base.ids import base62_12
from .collection import Collection
from .project import Person, Project
from .record import Space
from .schema import Schema
from .ulabel import ULabel


class DataMixin(models.Model):
    space: Space = ForeignKey(Space, PROTECT, default=1, db_default=1)
    feature = ForeignKey(
        Feature, null=True, blank=True, on_delete=models.CASCADE, related_name="+"
    )
    param = ForeignKey(
        Param, null=True, blank=True, on_delete=models.CASCADE, related_name="+"
    )
    row = IntegerField(help_text="Use -1 for result data")

    # Value fields
    value_int = models.BigIntegerField(null=True, blank=True)
    value_float = models.FloatField(null=True, blank=True)
    value_str = models.TextField(null=True, blank=True)
    value_datetime = models.DateTimeField(null=True, blank=True)
    value_ulabel = models.ForeignKey(
        ULabel, null=True, blank=True, on_delete=models.CASCADE, related_name="+"
    )
    value_person = models.ForeignKey(
        Person, null=True, blank=True, on_delete=models.CASCADE, related_name="+"
    )
    value_artifact = models.ForeignKey(
        Artifact, null=True, blank=True, on_delete=models.CASCADE, related_name="+"
    )
    value_collection = models.ForeignKey(
        Collection, null=True, blank=True, on_delete=models.CASCADE, related_name="+"
    )
    value_project = models.ForeignKey(
        Project, null=True, blank=True, on_delete=models.CASCADE, related_name="+"
    )
    value_json = models.JSONField(null=True, blank=True)

    class Meta:
        abstract = True

    def clean(self):
        # Validate feature/param mutual exclusivity
        if (self.feature is not None) == (self.param is not None):
            raise ValidationError("Exactly one of feature or param must be set")

        # Validate value fields
        values = [
            self.value_int,
            self.value_float,
            self.value_str,
            self.value_datetime,
            self.value_ulabel,
            self.value_artifact,
            self.value_json,
        ]
        non_null_count = sum(1 for v in values if v is not None)

        if non_null_count != 1:
            raise ValidationError("Exactly one value field must be set")


class RunData(BasicRecord, DataMixin):
    run = models.ForeignKey("Run", on_delete=models.CASCADE, related_name="_rundata")

    class Meta:
        constraints = [
            models.CheckConstraint(
                condition=(
                    models.Q(feature__isnull=False, param__isnull=True)
                    | models.Q(feature__isnull=True, param__isnull=False)
                ),
                name="run_data_feature_param_mutex",
            ),
            models.UniqueConstraint(
                fields=["run", "row", "feature", "param"], name="run_data_unique"
            ),
        ]
        indexes = [
            models.Index(fields=["run", "row"]),
            models.Index(fields=["feature"]),
            models.Index(fields=["param"]),
        ]


class FlexTable(Record, TracksRun, TracksUpdates):
    uid: str = CharField(
        editable=False, unique=True, max_length=12, db_index=True, default=base62_12
    )
    name = CharField()
    schema: Schema | None = ForeignKey(
        Schema, null=True, on_delete=models.SET_NULL, related_name="_tidytables"
    )
    type: Optional["FlexTable"] = ForeignKey(
        "self", PROTECT, null=True, related_name="records"
    )
    """Type of tidy table, e.g., `Cell`, `SampleSheet`, etc."""
    records: ULabel
    """Records of this type."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    description: str = CharField(null=True, db_index=True)
    """A description."""
    projects: Project = ManyToManyField(Project, related_name="_tidytables")
    ulabels: Project = ManyToManyField(ULabel, related_name="_tidytables")

    class Meta:
        indexes = [models.Index(fields=["uid"]), models.Index(fields=["name"])]


class FlexTableData(BasicRecord, DataMixin):
    tidytable = models.ForeignKey(
        FlexTable, on_delete=models.CASCADE, related_name="data"
    )

    class Meta:
        constraints = [
            models.CheckConstraint(
                condition=(
                    models.Q(feature__isnull=False, param__isnull=True)
                    | models.Q(feature__isnull=True, param__isnull=False)
                ),
                name="tidy_table_data_feature_param_mutex",
            ),
            models.UniqueConstraint(
                fields=["tidytable", "row", "feature", "param"],
                name="tidy_table_data_unique",
            ),
        ]
        indexes = [
            models.Index(fields=["tidytable", "row"]),
            models.Index(fields=["feature"]),
            models.Index(fields=["param"]),
        ]
