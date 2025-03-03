from datetime import date as DateType
from typing import Optional

from django.core.validators import RegexValidator
from django.db import models
from django.db.models import CASCADE, PROTECT

from lamindb.base.fields import (
    BigIntegerField,
    BooleanField,
    CharField,
    DateField,
    EmailField,
    ForeignKey,
    TextField,
    URLField,
)

from ..base.ids import base62_8, base62_12
from .artifact import Artifact
from .can_curate import CanCurate
from .collection import Collection
from .feature import Feature
from .record import BasicRecord, LinkORM, Record, ValidateFields
from .run import TracksRun, TracksUpdates
from .schema import Schema
from .transform import Transform
from .ulabel import ULabel


class Person(Record, CanCurate, TracksRun, TracksUpdates, ValidateFields):
    """Persons.

    This registry is distinct from `User` and purely exists for project management.

    You'll soon be able to conveniently create persons from users.

    Example:
        >>> person = Person(
        ...     name="Jane Doe",
        ...     email="jane.doe@example.com",
        ...     internal=True,
        ... ).save()
    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, max_length=8, db_index=True, default=base62_8
    )
    """Universal id, valid across DB instances."""
    name: str = CharField(db_index=True)
    """Name of the person (forename(s) lastname)."""
    email: str | None = EmailField(null=True, default=None)
    """Email of the person."""
    external: bool = BooleanField(default=True, db_index=True)
    """Whether the person is external to the organization."""


class Reference(Record, CanCurate, TracksRun, TracksUpdates, ValidateFields):
    """References such as internal studies, papers, documents, or URLs.

    Example:
        >>> reference = Reference(
        ...     name="A Paper Title",
        ...     abbr="APT",
        ...     url="https://doi.org/10.1000/xyz123",
        ...     pubmed_id=12345678,
        ...     doi="10.1000/xyz123",
        ...     description="Good paper.",
        ...     text="Some text I want to be searchable.",
        ...     date=date(2023, 11, 21),
        ... ).save()
    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, max_length=12, db_index=True, default=base62_12
    )
    """Universal id, valid across DB instances."""
    name: str = CharField(db_index=True)
    """Title or name of the reference document."""
    abbr: str | None = CharField(
        max_length=32,
        db_index=True,
        null=True,
    )
    """An abbreviation for the reference."""
    type: Optional["Reference"] = ForeignKey(
        "self", PROTECT, null=True, related_name="records"
    )
    """Type of reference (e.g., 'Study', 'Paper', 'Preprint').

    Allows to group reference by type, e.g., internal studies vs. all papers etc.
    """
    records: "Reference"
    """Records of this type."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    url: str | None = URLField(null=True)
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
    description: str | None = CharField(null=True, db_index=True)
    """Description of the reference."""
    text: str | None = TextField(null=True)
    """Abstract or full text of the reference to make it searchable."""
    date: DateType | None = DateField(null=True, default=None)
    """Date of creation or publication of the reference."""
    authors: Person = models.ManyToManyField(Person, related_name="references")
    """All people associated with this reference."""
    artifacts: Artifact = models.ManyToManyField(
        Artifact, through="ArtifactReference", related_name="references"
    )
    """Artifacts associated with this reference."""
    transforms: Artifact = models.ManyToManyField(
        Transform, through="TransformReference", related_name="references"
    )
    """Transforms associated with this reference."""
    collections: Artifact = models.ManyToManyField(
        Collection, through="CollectionReference", related_name="references"
    )
    """Collections associated with this reference."""


class Project(Record, CanCurate, TracksRun, TracksUpdates, ValidateFields):
    """Projects.

    Example:
        >>> project = Project(
        ...     name="My Project Name",
        ...     abbr="MPN",
        ...     url="https://example.com/my_project",
        ... ).save()
    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, max_length=12, db_index=True, default=base62_12
    )
    """Universal id, valid across DB instances."""
    name: str = CharField(db_index=True)
    """Title or name of the Project."""
    type: Optional["Project"] = ForeignKey(
        "self", PROTECT, null=True, related_name="records"
    )
    """Type of project (e.g., 'Program', 'Project', 'GithubIssue', 'Task')."""
    records: "Project"
    """Records of this type."""
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
    parents: "Project" = models.ManyToManyField(
        "self", symmetrical=False, related_name="children"
    )
    """Parent projects, the super-projects owning this project."""
    children: "Project"
    """Child projects, the sub-projects owned by this project.

    Reverse accessor for `.parents`.
    """
    predecessors: "Project" = models.ManyToManyField(
        "self", symmetrical=False, related_name="successors"
    )
    """The preceding projects required by this project."""
    successors: "Project"
    """The succeeding projects requiring this project.

    Reverse accessor for `.predecessors`.
    """
    people: Person = models.ManyToManyField(
        Person, through="PersonProject", related_name="projects"
    )
    """People associated with this project."""
    artifacts: Artifact = models.ManyToManyField(
        Artifact, through="ArtifactProject", related_name="projects"
    )
    """Artifacts associated with this Project."""
    transforms: Transform = models.ManyToManyField(
        Transform, through="TransformProject", related_name="projects"
    )
    """Transforms associated with this project."""
    ulabels: ULabel = models.ManyToManyField(
        ULabel, through="ULabelProject", related_name="projects"
    )
    """Transforms associated with this project."""
    features: ULabel = models.ManyToManyField(
        Feature, through="FeatureProject", related_name="projects"
    )
    """Transforms associated with this project."""
    schemas: ULabel = models.ManyToManyField(
        Schema, through="SchemaProject", related_name="projects"
    )
    """Schemas associated with this project."""
    collections: Collection = models.ManyToManyField(
        Collection, through="CollectionProject", related_name="projects"
    )
    """Collections associated with this project."""
    references: "Reference" = models.ManyToManyField(
        "Reference", related_name="projects"
    )
    """References associated with this project."""
    _status_code: int = models.SmallIntegerField(default=0, db_index=True)
    """Status code."""


class ArtifactProject(BasicRecord, LinkORM, TracksRun):
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
        # can have the same label linked to the same artifact if the feature is different
        unique_together = ("artifact", "project", "feature")


class TransformProject(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    transform: Transform = ForeignKey(Transform, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_transform")

    class Meta:
        unique_together = ("transform", "project")


class CollectionProject(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="links_project"
    )
    project: Project = ForeignKey(Project, PROTECT, related_name="links_collection")

    class Meta:
        unique_together = ("collection", "project")


class ULabelProject(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    ulabel: Transform = ForeignKey(ULabel, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_ulabel")

    class Meta:
        unique_together = ("ulabel", "project")


class PersonProject(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    person: Transform = ForeignKey(Person, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_person")
    role: str | None = CharField(null=True, default=None)

    class Meta:
        unique_together = ("person", "project")


class FeatureProject(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    feature: Feature = ForeignKey(Feature, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_feature")

    class Meta:
        unique_together = ("feature", "project")


class SchemaProject(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    schema: Schema = ForeignKey(Schema, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_schema")

    class Meta:
        unique_together = ("schema", "project")


class ArtifactReference(BasicRecord, LinkORM, TracksRun):
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
        # can have the same label linked to the same artifact if the feature is different
        unique_together = ("artifact", "reference", "feature")


class TransformReference(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    transform: Transform = ForeignKey(
        Transform, CASCADE, related_name="links_reference"
    )
    reference: Reference = ForeignKey(
        Reference, PROTECT, related_name="links_transform"
    )

    class Meta:
        unique_together = ("transform", "reference")


class CollectionReference(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="links_reference"
    )
    reference: Reference = ForeignKey(
        Reference, PROTECT, related_name="links_collection"
    )

    class Meta:
        unique_together = ("collection", "reference")
