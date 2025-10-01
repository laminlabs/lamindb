from datetime import datetime

from django.db import models
from django.db.models import CASCADE, CharField, DateTimeField, ForeignKey, TextField

from .artifact import Artifact
from .run import User
from .sqlrecord import BaseSQLRecord, IsVersioned


class RootPage(BaseSQLRecord, IsVersioned):
    """A root page for every registry that can appear at the top of the registry root page in the GUI."""

    _len_full_uid: int = 20
    _len_stem_uid: int = 16

    class Meta:
        app_label = "lamindb"

    id = models.BigAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=_len_full_uid
    )
    """Universal id."""
    table_name: str = CharField(max_length=255, db_index=True)
    """The SQL table name."""
    content: str = TextField()
    """Markdown content of the page."""
    hash: str | None = CharField(max_length=22, db_index=True, null=True)
    """Content hash of the page."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=None, related_name="+", null=True
    )
    """Creator of page."""


class ArtifactPage(BaseSQLRecord, IsVersioned):
    """An unstructured notes page that can be attached to an artifact."""

    _len_full_uid: int = 20
    _len_stem_uid: int = 16

    class Meta:
        app_label = "lamindb"

    id = models.BigAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=_len_full_uid
    )
    """Universal id."""
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="pages", null=True)
    """The artifact to which the page is attached."""
    content: str = TextField()
    """Markdown content of the page."""
    hash: str | None = CharField(max_length=22, db_index=True, null=True)
    """Content hash of the page."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=None, related_name="+", null=True
    )
    """Creator of page."""
