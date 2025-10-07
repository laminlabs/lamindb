from datetime import datetime

from django.db import models
from django.db.models import CASCADE, CharField, DateTimeField, ForeignKey, TextField

from .artifact import Artifact
from .collection import Collection
from .feature import Feature
from .project import Project
from .record import Record
from .run import Run, User
from .schema import Schema
from .sqlrecord import BaseSQLRecord, Branch, IsVersioned, Space
from .transform import Transform


class BlockMixin(IsVersioned):
    class Meta:
        abstract = True

    _len_full_uid: int = 20
    _len_stem_uid: int = 16

    id = models.BigAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=_len_full_uid
    )
    """Universal id."""
    content: str = TextField()
    """Content of the block."""
    hash: str = CharField(max_length=22, db_index=True, null=True)
    """Content hash of the block."""
    type: str = CharField(
        max_length=22, db_index=True, default="mdpage", db_default="mdpage"
    )
    """Type of the block.

    Only current option is: "mdpage" (markdown page).
    """
    vertical_pos: float = models.FloatField(default=0, db_default=0, db_index=True)
    """The vertical position of the block in the GUI as a float.

    It can be used to act as a percentage to define absolute positioning or for ordering.
    """
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=None, related_name="+", null=True
    )
    """Creator of block."""


class RootBlock(BlockMixin, BaseSQLRecord):
    """A root block for every registry that can appear at the top of the registry root block in the GUI."""

    class Meta:
        app_label = "lamindb"

    name: str = CharField(max_length=255, db_index=True)
    """The entity for which we want to create a block.

    Conventions are mostly to take the SQL table name:
        name = "instance" means instance
        name = "lamindb_artifact" means artifact
        name = "lamindb_transform" means transform
        name = "bionty_celltype" means bionty cell type
    """


class RecordBlock(BlockMixin, BaseSQLRecord):
    """An unstructured notes block that can be attached to a record."""

    class Meta:
        app_label = "lamindb"

    record: Record = ForeignKey(Record, CASCADE, related_name="blocks")
    """The record to which the block is attached."""


class ArtifactBlock(BlockMixin, BaseSQLRecord):
    """An unstructured notes block that can be attached to an artifact."""

    class Meta:
        app_label = "lamindb"

    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="blocks")
    """The artifact to which the block is attached."""


class TransformBlock(BlockMixin, BaseSQLRecord):
    """An unstructured notes block that can be attached to a transform."""

    class Meta:
        app_label = "lamindb"

    transform: Transform = ForeignKey(
        Transform, CASCADE, related_name="blocks", null=True
    )
    """The transform to which the block is attached."""
    line_number: int | None = models.IntegerField(null=True)
    """The line number in the source code to which the block belongs."""


class RunBlock(BlockMixin, BaseSQLRecord):
    """An unstructured notes block that can be attached to a run."""

    class Meta:
        app_label = "lamindb"

    run: Run = ForeignKey(Run, CASCADE, related_name="blocks")
    """The run to which the block is attached."""


class CollectionBlock(BlockMixin, BaseSQLRecord):
    """An unstructured notes block that can be attached to a collection."""

    class Meta:
        app_label = "lamindb"

    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="blocks", null=True
    )
    """The collection to which the block is attached."""


class SchemaBlock(BlockMixin, BaseSQLRecord):
    """An unstructured notes block that can be attached to a schema."""

    class Meta:
        app_label = "lamindb"

    schema: Schema = ForeignKey(Schema, CASCADE, related_name="blocks")
    """The schema to which the block is attached."""


class FeatureBlock(BlockMixin, BaseSQLRecord):
    """An unstructured notes block that can be attached to a feature."""

    class Meta:
        app_label = "lamindb"

    feature: Feature = ForeignKey(Feature, CASCADE, related_name="blocks")
    """The feature to which the block is attached."""


class ProjectBlock(BlockMixin, BaseSQLRecord):
    """An unstructured notes block that can be attached to a project."""

    class Meta:
        app_label = "lamindb"

    project: Project = ForeignKey(Project, CASCADE, related_name="blocks")
    """The project to which the block is attached."""


class SpaceBlock(BlockMixin, BaseSQLRecord):
    """An unstructured notes block that can be attached to a space."""

    class Meta:
        app_label = "lamindb"

    space: Space = ForeignKey(Space, CASCADE, related_name="blocks")
    """The space to which the block is attached."""


class BranchBlock(BlockMixin, BaseSQLRecord):
    """An unstructured notes block that can be attached to a branch."""

    class Meta:
        app_label = "lamindb"

    branch: Branch = ForeignKey(Branch, CASCADE, related_name="blocks")
    """The branch to which the block is attached."""
