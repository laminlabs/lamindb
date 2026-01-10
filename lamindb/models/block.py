from datetime import datetime
from typing import Any

from django.db import models
from django.db.models import (
    CASCADE,
    CharField,
    DateTimeField,
    ForeignKey,
    JSONField,
    TextField,
)

from ..base.uids import base62_16
from .artifact import Artifact
from .collection import Collection
from .feature import Feature
from .project import Project
from .record import Record
from .run import Run, User
from .schema import Schema
from .sqlrecord import BaseSQLRecord, Branch, IsVersioned, Space, SQLRecord
from .transform import Transform


class BaseBlock(IsVersioned):
    class Meta:
        abstract = True

    _len_full_uid: int = 20
    _len_stem_uid: int = 16

    id = models.BigAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False,
        unique=True,
        db_index=True,
        max_length=_len_full_uid,
        default=base62_16,
    )
    """Universal id."""
    content: str = TextField()
    """Content of the block."""
    hash: str = CharField(max_length=22, db_index=True, null=True)
    """Content hash of the block."""
    kind: str = CharField(
        max_length=22, db_index=True, default="mdpage", db_default="mdpage"
    )
    """Kind of block.

    Only current option is: "mdpage" (markdown page).
    """
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=None, related_name="+", null=True
    )
    """Creator of block."""
    _aux: dict[str, Any] | None = JSONField(default=None, db_default=None, null=True)
    """Auxiliary field for dictionary-like metadata."""


class Block(BaseBlock, SQLRecord):
    """A root block for every registry that can appear at the top of the registry root block in the GUI."""

    class Meta:
        app_label = "lamindb"

    # same key as in transform/artifact/collection
    key: str = CharField(max_length=1024, db_index=True)
    """The key for which we want to create a block."""
    projects: Project
    """Projects that annotate this block."""


class RecordBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a record."""

    class Meta:
        app_label = "lamindb"

    record: Record = ForeignKey(Record, CASCADE, related_name="ablocks")
    """The record to which the block is attached."""


class ArtifactBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to an artifact."""

    class Meta:
        app_label = "lamindb"

    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="ablocks")
    """The artifact to which the block is attached."""


class TransformBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a transform."""

    class Meta:
        app_label = "lamindb"

    transform: Transform = ForeignKey(
        Transform, CASCADE, related_name="ablocks", null=True
    )
    """The transform to which the block is attached."""
    line_number: int | None = models.IntegerField(null=True)
    """The line number in the source code to which the block belongs."""


class RunBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a run."""

    class Meta:
        app_label = "lamindb"

    run: Run = ForeignKey(Run, CASCADE, related_name="ablocks")
    """The run to which the block is attached."""


class CollectionBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a collection."""

    class Meta:
        app_label = "lamindb"

    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="ablocks", null=True
    )
    """The collection to which the block is attached."""


class SchemaBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a schema."""

    class Meta:
        app_label = "lamindb"

    schema: Schema = ForeignKey(Schema, CASCADE, related_name="ablocks")
    """The schema to which the block is attached."""


class FeatureBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a feature."""

    class Meta:
        app_label = "lamindb"

    feature: Feature = ForeignKey(Feature, CASCADE, related_name="ablocks")
    """The feature to which the block is attached."""


class ProjectBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a project."""

    class Meta:
        app_label = "lamindb"

    project: Project = ForeignKey(Project, CASCADE, related_name="ablocks")
    """The project to which the block is attached."""


class SpaceBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a space."""

    class Meta:
        app_label = "lamindb"

    space: Space = ForeignKey(Space, CASCADE, related_name="ablocks")
    """The space to which the block is attached."""


class BranchBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a branch."""

    class Meta:
        app_label = "lamindb"

    branch: Branch = ForeignKey(Branch, CASCADE, related_name="ablocks")
    """The branch to which the block is attached."""


class ULabelBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a ulabel."""

    class Meta:
        app_label = "lamindb"

    ulabel = ForeignKey("ULabel", CASCADE, related_name="ablocks")
    """The ulabel to which the block is attached."""
