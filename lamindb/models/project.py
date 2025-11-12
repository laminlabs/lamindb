from __future__ import annotations

from typing import TYPE_CHECKING, overload

from django.core.validators import RegexValidator
from django.db import models
from django.db.models import CASCADE, PROTECT

from lamindb.base.fields import (
    BigIntegerField,
    BooleanField,
    CharField,
    DateField,
    DateTimeField,
    ForeignKey,
    TextField,
    URLField,
)
from lamindb.base.users import current_user_id

from ..base.ids import base62_12
from .artifact import Artifact
from .can_curate import CanCurate
from .collection import Collection
from .feature import Feature
from .record import Record
from .run import Run, TracksRun, TracksUpdates, User
from .schema import Schema
from .sqlrecord import BaseSQLRecord, IsLink, SQLRecord, ValidateFields
from .transform import Transform
from .ulabel import ULabel

if TYPE_CHECKING:
    from datetime import date as DateType
    from datetime import datetime

    from .block import ProjectBlock


class Reference(SQLRecord, CanCurate, TracksRun, TracksUpdates, ValidateFields):
    """References such as internal studies, papers, documents, or URLs.

    Example:

        ::
            reference = Reference(
                name="A Paper Title",
                abbr="APT",
                url="https://doi.org/10.1000/xyz123",
                pubmed_id=12345678,
                doi="10.1000/xyz123",
                description="Good paper.",
                text="Some text I want to be searchable.",
                date=date(2023, 11, 21),
            ).save()
    """

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False
        app_label = "lamindb"
        constraints = [
            models.UniqueConstraint(
                fields=["name", "type", "space"],
                name="unique_reference_name_type_space",
                condition=~models.Q(branch_id=-1),
            )
            # also see raw SQL constraints for `is_type` and `type` FK validity in migrations
        ]

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, max_length=12, db_index=True, default=base62_12
    )
    """Universal id, valid across DB instances."""
    name: str = CharField(db_index=True)
    """Title or name of the reference document."""
    description: str | None = TextField(null=True)
    """A description."""
    type: Reference | None = ForeignKey(
        "self", PROTECT, null=True, related_name="references"
    )
    """Type of reference (e.g., 'Study', 'Paper', 'Preprint').

    Allows to group reference by type, e.g., internal studies vs. all papers etc.
    """
    references: Reference
    """References of this type (can only be non-empty if `is_type` is `True`)."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    abbr: str | None = CharField(
        max_length=32,
        db_index=True,
        null=True,
    )
    """An abbreviation for the reference."""
    url: str | None = URLField(null=True, db_index=True)
    """URL linking to the reference."""
    pubmed_id: int | None = BigIntegerField(null=True, db_index=True)
    """A PudMmed ID."""
    doi: str | None = CharField(
        null=True,
        db_index=True,
        validators=[
            RegexValidator(
                regex=r"^(?:https?://(?:dx\.)?doi\.org/|doi:|DOI:)?10\.\d+/.*$",
                message="Must be a DOI (e.g., 10.1000/xyz123 or https://doi.org/10.1000/xyz123)",
            )
        ],
    )
    """Digital Object Identifier (DOI) for the reference."""
    text: str | None = TextField(null=True)
    """Abstract or full text of the reference to make it searchable."""
    date: DateType | None = DateField(null=True, default=None)
    """Date of creation or publication of the reference."""
    artifacts: Artifact = models.ManyToManyField(
        Artifact, through="ArtifactReference", related_name="references"
    )
    """Annotated artifacts."""
    transforms: Transform = models.ManyToManyField(
        Transform, through="TransformReference", related_name="references"
    )
    """Annotated transforms."""
    collections: Collection = models.ManyToManyField(
        Collection, through="CollectionReference", related_name="references"
    )
    """Annotated collections."""
    linked_in_records: Record = models.ManyToManyField(
        Record, through="RecordReference", related_name="linked_references"
    )
    """Linked in records."""
    records: Record = models.ManyToManyField(
        Record, through="ReferenceRecord", related_name="references"
    )
    """Annotated records."""
    projects: Project
    """Projects that annotate this reference."""

    @overload
    def __init__(
        self,
        name: str,
        type: Reference | None = None,
        is_type: bool = False,
        abbr: str | None = None,
        url: str | None = None,
        pubmed_id: int | None = None,
        doi: str | None = None,
        description: str | None = None,
        text: str | None = None,
        date: DateType | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class Project(SQLRecord, CanCurate, TracksRun, TracksUpdates, ValidateFields):
    """Projects to label artifacts, transforms, records, and runs.

    Example:

        ::

            project = Project(
            name="My Project Name",
            abbr="MPN",
            url="https://example.com/my_project",
            ).save()
            artifact.projects.add(project)  # <-- labels the artifact with the project
            ln.track(project=project)       # <-- automtically labels entities during the run

    """

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False
        app_label = "lamindb"
        constraints = [
            models.UniqueConstraint(
                fields=["name", "type", "space"],
                name="unique_project_name_type_space",
                condition=~models.Q(branch_id=-1),
            )
            # also see raw SQL constraints for `is_type` and `type` FK validity in migrations
        ]

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, max_length=12, db_index=True, default=base62_12
    )
    """Universal id, valid across DB instances."""
    name: str = CharField(db_index=True)
    """Title or name of the Project."""
    description: str | None = TextField(null=True)
    """A description."""
    type: Project | None = ForeignKey(
        "self", PROTECT, null=True, related_name="projects"
    )
    """Type of project (e.g., 'Program', 'Project', 'GithubIssue', 'Task')."""
    projects: Project
    """Projects of this type (can only be non-empty if `is_type` is `True`)."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    abbr: str | None = CharField(max_length=32, db_index=True, null=True)
    """An abbreviation."""
    url: str | None = URLField(max_length=255, null=True, default=None)
    """A URL."""
    start_date: DateType | None = DateField(null=True, default=None)
    """Date of start of the project."""
    end_date: DateType | None = DateField(null=True, default=None)
    """Date of start of the project."""
    parents: Project = models.ManyToManyField(
        "self", symmetrical=False, related_name="children"
    )
    """Parent projects, the super-projects owning this project."""
    children: Project
    """Child projects, the sub-projects owned by this project.

    Reverse accessor for `.parents`.
    """
    predecessors: Project = models.ManyToManyField(
        "self", symmetrical=False, related_name="successors"
    )
    """The preceding projects required by this project."""
    successors: Project
    """The succeeding projects requiring this project.

    Reverse accessor for `.predecessors`.
    """
    artifacts: Artifact = models.ManyToManyField(
        Artifact, through="ArtifactProject", related_name="projects"
    )
    """Annotated artifacts."""
    transforms: Transform = models.ManyToManyField(
        Transform, through="TransformProject", related_name="projects"
    )
    """Annotated transforms."""
    runs: Run = models.ManyToManyField(
        Run, through="RunProject", related_name="projects"
    )
    """Annotated runs."""
    ulabels: ULabel = models.ManyToManyField(
        ULabel, through="ULabelProject", related_name="projects"
    )
    """Annotated ulabels."""
    features: Feature = models.ManyToManyField(
        Feature, through="FeatureProject", related_name="projects"
    )
    """Annotated features."""
    schemas: Schema = models.ManyToManyField(
        Schema, through="SchemaProject", related_name="projects"
    )
    """Annotated schemas."""
    linked_in_records: Record = models.ManyToManyField(
        Record, through="RecordProject", related_name="linked_projects"
    )
    """Linked in records."""
    records: Record = models.ManyToManyField(
        Record, through="ProjectRecord", related_name="projects"
    )
    """Annotated records."""
    collections: Collection = models.ManyToManyField(
        Collection, through="CollectionProject", related_name="projects"
    )
    """Annotated collections."""
    references: Reference = models.ManyToManyField("Reference", related_name="projects")
    """Annotated references."""
    _status_code: int = models.SmallIntegerField(default=0, db_index=True)
    """Status code."""
    blocks: ProjectBlock
    """Blocks that annotate this project."""

    @overload
    def __init__(
        self,
        name: str,
        type: Project | None = None,
        is_type: bool = False,
        abbr: str | None = None,
        url: str | None = None,
        start_date: DateType | None = None,
        end_date: DateType | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class ArtifactProject(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_artifact")
    feature: Feature | None = ForeignKey(
        Feature,
        PROTECT,
        null=True,
        default=None,
        related_name="links_artifactproject",
    )
    label_ref_is_name: bool | None = BooleanField(null=True, default=None)
    feature_ref_is_name: bool | None = BooleanField(null=True, default=None)

    class Meta:
        app_label = "lamindb"
        # can have the same label linked to the same artifact if the feature is different
        unique_together = ("artifact", "project", "feature")


class RunProject(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    run: Run = ForeignKey(Run, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_run")
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

    class Meta:
        app_label = "lamindb"
        unique_together = ("run", "project")


class TransformProject(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    transform: Transform = ForeignKey(Transform, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_transform")

    class Meta:
        app_label = "lamindb"
        unique_together = ("transform", "project")


class CollectionProject(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="links_project"
    )
    project: Project = ForeignKey(Project, PROTECT, related_name="links_collection")

    class Meta:
        app_label = "lamindb"
        unique_together = ("collection", "project")


class ULabelProject(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    ulabel: ULabel = ForeignKey(ULabel, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_ulabel")

    class Meta:
        app_label = "lamindb"
        unique_together = ("ulabel", "project")


class FeatureProject(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    feature: Feature = ForeignKey(Feature, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_feature")

    class Meta:
        app_label = "lamindb"
        unique_together = ("feature", "project")


class SchemaProject(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    schema: Schema = ForeignKey(Schema, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_schema")

    class Meta:
        app_label = "lamindb"
        unique_together = ("schema", "project")


# for annotation of records with references, RecordReference is for storing reference values
class ReferenceRecord(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    reference: Reference = ForeignKey(Reference, PROTECT, related_name="links_record")
    feature: Feature | None = ForeignKey(
        Feature,
        PROTECT,
        null=True,
        default=None,
        related_name="links_referencerecord",
    )
    record: Record = ForeignKey(Record, CASCADE, related_name="links_reference")

    class Meta:
        app_label = "lamindb"
        unique_together = ("reference", "feature", "record")


class RecordReference(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_reference")
    feature: Feature = ForeignKey(
        Feature, PROTECT, related_name="links_recordreference"
    )
    value: Reference = ForeignKey(Reference, PROTECT, related_name="links_in_record")

    class Meta:
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for annotation of records with projects, RecordProject is for storing project values
class ProjectRecord(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    project: Project = ForeignKey(Project, PROTECT, related_name="links_record")
    feature: Feature | None = ForeignKey(
        Feature,
        PROTECT,
        null=True,
        default=None,
        related_name="links_projectrecord",
    )
    record: Record = ForeignKey(Record, CASCADE, related_name="links_project")

    class Meta:
        app_label = "lamindb"
        unique_together = ("project", "feature", "record")


class RecordProject(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_project")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordproject")
    value: Project = ForeignKey(Project, PROTECT, related_name="links_in_record")

    class Meta:
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


class ArtifactReference(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="links_reference")
    reference: Reference = ForeignKey(Reference, PROTECT, related_name="links_artifact")
    feature: Feature | None = ForeignKey(
        Feature,
        PROTECT,
        null=True,
        default=None,
        related_name="links_artifactreference",
    )
    label_ref_is_name: bool | None = BooleanField(null=True, default=None)
    feature_ref_is_name: bool | None = BooleanField(null=True, default=None)

    class Meta:
        app_label = "lamindb"
        unique_together = ("artifact", "reference", "feature")


class TransformReference(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    transform: Transform = ForeignKey(
        Transform, CASCADE, related_name="links_reference"
    )
    reference: Reference = ForeignKey(
        Reference, PROTECT, related_name="links_transform"
    )

    class Meta:
        app_label = "lamindb"
        unique_together = ("transform", "reference")


class CollectionReference(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="links_reference"
    )
    reference: Reference = ForeignKey(
        Reference, PROTECT, related_name="links_collection"
    )

    class Meta:
        app_label = "lamindb"
        unique_together = ("collection", "reference")
