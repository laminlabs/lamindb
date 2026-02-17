from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, overload

from django.db import models
from django.db.models import (
    CASCADE,
    PROTECT,
    CharField,
    DateTimeField,
    ForeignKey,
    JSONField,
    Q,
    TextField,
)
from lamin_utils import logger
from lamindb_setup.core.hashing import hash_string

from ..base.uids import base62_16
from ._is_versioned import create_uid, process_revises
from .artifact import Artifact
from .collection import Collection
from .feature import Feature
from .project import Project
from .record import Record
from .run import Run, User
from .schema import Schema
from .sqlrecord import (
    BaseSQLRecord,
    Branch,
    IsVersioned,
    Space,
    SQLRecord,
    init_self_from_db,
    update_attributes,
)
from .transform import Transform

if TYPE_CHECKING:
    from datetime import datetime

    from .query_manager import RelatedManager

_VERSIONED_ATTACHED_KINDS = ("readme",)  # only readme is versioned; comment is not
_VALID_BLOCK_KINDS: tuple[str, ...] = ("readme", "comment")


def _init_versioned_attached_block(
    self: BaseBlock,
    fk_field_name: str,
    *args: Any,
    allowed_extra: tuple[str, ...] = (),
    **kwargs: Any,
) -> None:
    cls = type(self)
    if len(args) == len(self._meta.concrete_fields):
        super(cls, self).__init__(*args, **kwargs)
        return None
    if args:
        raise ValueError(
            f"Please only use keyword arguments to construct a {cls.__name__}"
        )
    fk_value = kwargs.pop(fk_field_name, None)
    content = kwargs.pop("content", None)
    kind = kwargs.pop("kind", None)
    version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
    revises = kwargs.pop("revises", None)
    using = kwargs.pop("using", None)
    uid = kwargs.pop("uid", None) if "uid" in kwargs else None
    extra_kwargs = {k: kwargs.pop(k) for k in allowed_extra if k in kwargs}
    allowed = {
        fk_field_name,
        "content",
        "kind",
        "version",
        "version_tag",
        "revises",
        "using",
        "uid",
        *allowed_extra,
    }
    if kwargs:
        raise ValueError(
            f"Only {', '.join(sorted(allowed))} can be passed, but you passed: {kwargs}"
        )
    if fk_value is None:
        raise ValueError(f"{fk_field_name} is required for {cls.__name__}")
    if kind is None:
        raise ValueError(
            f"kind is required for {cls.__name__}; use 'readme' or 'comment'"
        )
    if kind not in _VALID_BLOCK_KINDS:
        raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")

    if kind == "comment":
        if revises is not None:
            raise ValueError(
                "revises is not allowed for kind='comment'; comments are not versioned"
            )
        new_uid, _ = create_uid(
            revises=None,
            version_tag=version_tag,
            n_full_id=cls._len_full_uid,
        )
        block_hash = hash_string(content) if content else None
        super(cls, self).__init__(
            uid=new_uid,
            content=content or "",
            hash=block_hash,
            kind=kind,
            version_tag=version_tag,
            revises=None,
            **{fk_field_name: fk_value},
            **extra_kwargs,
        )
        return None
    # kind == "readme" (versioned)
    if revises is None and fk_value is not None:
        candidate_for_revises = (
            cls.objects.using(using)
            .filter(
                **{fk_field_name: fk_value},
                kind=kind,
                is_latest=True,
            )
            .order_by("-created_at")
            .first()
        )
        if candidate_for_revises is not None:
            revises = candidate_for_revises
            content_blank = getattr(revises, "content", None) in (None, "")
            if content_blank:
                logger.important(
                    "no content was yet saved, returning existing "
                    f"block with same {fk_field_name} and kind"
                )
                uid = revises.uid
    if revises is not None and uid is not None and uid == revises.uid:
        init_self_from_db(self, revises)
        update_attributes(self, {})
        return None
    new_uid, revises = create_uid(
        revises=revises,
        version_tag=version_tag,
        n_full_id=cls._len_full_uid,
    )
    if uid is None:
        uid = new_uid
    block_hash = hash_string(content) if content else None
    super(cls, self).__init__(
        uid=uid,
        content=content or "",
        hash=block_hash,
        kind=kind,
        version_tag=version_tag,
        revises=revises,
        **{fk_field_name: fk_value},
        **extra_kwargs,
    )


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
        max_length=22, db_index=True, default="readme", db_default="readme"
    )
    """The kind of block: "readme" (readme-type markdown page) or "comment"."""
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
    """An experimental markdown block for anything: issues, standalone markdown pages, comments, etc.

    The `Block` model is experimental and may change in the future.
    """

    class Meta:
        app_label = "lamindb"

    # same key as in transform/artifact/collection
    key: str | None = CharField(max_length=1024, db_index=True, null=True)
    """The key for which we want to create a block."""
    anchor: Block | None = ForeignKey(
        "Block", PROTECT, related_name="children", null=True
    )
    """The anchor of this block.

    For a comment, could be the issue on which the comment is attached.

    For a sub-post, could be the parent post.
    """
    projects: RelatedManager[Project]
    """Projects that annotate this block."""
    anchors: RelatedManager[Block]
    """This block anchors these blocks."""

    @overload
    def __init__(
        self,
        key: str | None = None,
        content: str | None = None,
        kind: Literal["readme"] = ...,
        version: str | None = None,
        revises: Block | None = None,
        anchor: Block | None = None,
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
        if args:
            raise ValueError("Please only use keyword arguments to construct a Block")
        key = kwargs.pop("key", None)
        content = kwargs.pop("content", None)
        revises = kwargs.pop("revises", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        kind = kwargs.pop("kind", None)
        anchor = kwargs.pop("anchor", None)
        using = kwargs.pop("using", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        branch = kwargs.pop("branch", None)
        branch_id = kwargs.pop("branch_id", 1)
        space = kwargs.pop("space", None)
        space_id = kwargs.pop("space_id", 1)
        if kwargs:
            raise ValueError(
                "Only key, content, kind, version, revises, anchor "
                f"can be passed, but you passed: {kwargs}"
            )
        if kind != "readme":
            raise ValueError("Only kind = 'readme' is supported for block.")
        assert key.startswith("__lamindb_..."), "key must start with '__lamindb_...'"
        if revises is None:
            if uid is not None:
                revises = (
                    Block.objects.using(using)
                    .filter(
                        uid__startswith=uid[:-4],
                        is_latest=True,
                    )
                    .order_by("-created_at")
                    .first()
                )
            elif key is not None:
                candidate_for_revises = (
                    Block.objects.using(using)
                    .filter(
                        ~Q(branch_id=-1),
                        key=key,
                        is_latest=True,
                    )
                    .order_by("-created_at")
                    .first()
                )
                if candidate_for_revises is not None:
                    revises = candidate_for_revises
                    content_blank = getattr(candidate_for_revises, "content", None) in (
                        None,
                        "",
                    )
                    if content_blank:
                        logger.important(
                            "no content was yet saved, returning existing "
                            "block with same key"
                        )
                        uid = revises.uid
        if revises is not None and uid is not None and uid == revises.uid:
            if revises.key != key:
                logger.warning("ignoring inconsistent key")
            init_self_from_db(self, revises)
            update_attributes(self, {})
            return None
        if revises is not None and key is not None and revises.key != key:
            logger.important(f"renaming block {revises.key} to {key}")
        new_uid, version_tag, key, _, revises = process_revises(
            revises, version_tag, key, None, Block
        )
        if uid is None:
            uid = new_uid
        block_hash = None
        if content is not None:
            block_hash = hash_string(content)
            block_candidate = Block.objects.filter(
                ~Q(branch_id=-1),
                hash=block_hash,
                is_latest=True,
            ).first()
            if block_candidate is not None:
                init_self_from_db(self, block_candidate)
                update_attributes(self, {})
                if key is not None and block_candidate.key != key:
                    logger.warning(
                        f"key {self.key} on existing block differs from "
                        f"passed key {key}, keeping original key"
                    )
                return None
        super().__init__(
            uid=uid,
            key=key,
            content=content or "",
            kind=kind,
            version_tag=version_tag,
            hash=block_hash,
            revises=revises,
            anchor=anchor,
            branch=branch,
            branch_id=branch_id,
            space=space,
            space_id=space_id,
        )


class RecordBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a record."""

    class Meta:
        app_label = "lamindb"

    record: Record = ForeignKey(Record, CASCADE, related_name="ablocks")
    """The record to which the block is attached."""

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(self, "record", *args, **kwargs)


class ArtifactBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to an artifact."""

    class Meta:
        app_label = "lamindb"

    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="ablocks")
    """The artifact to which the block is attached."""

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(self, "artifact", *args, **kwargs)


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

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(
            self, "transform", *args, allowed_extra=("line_number",), **kwargs
        )


class RunBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a run."""

    class Meta:
        app_label = "lamindb"

    run: Run = ForeignKey(Run, CASCADE, related_name="ablocks")
    """The run to which the block is attached."""

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(self, "run", *args, **kwargs)


class CollectionBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a collection."""

    class Meta:
        app_label = "lamindb"

    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="ablocks", null=True
    )
    """The collection to which the block is attached."""

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(self, "collection", *args, **kwargs)


class SchemaBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a schema."""

    class Meta:
        app_label = "lamindb"

    schema: Schema = ForeignKey(Schema, CASCADE, related_name="ablocks")
    """The schema to which the block is attached."""

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(self, "schema", *args, **kwargs)


class FeatureBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a feature."""

    class Meta:
        app_label = "lamindb"

    feature: Feature = ForeignKey(Feature, CASCADE, related_name="ablocks")
    """The feature to which the block is attached."""

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(self, "feature", *args, **kwargs)


class ProjectBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a project."""

    class Meta:
        app_label = "lamindb"

    project: Project = ForeignKey(Project, CASCADE, related_name="ablocks")
    """The project to which the block is attached."""

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(self, "project", *args, **kwargs)


class SpaceBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a space."""

    class Meta:
        app_label = "lamindb"

    space: Space = ForeignKey(Space, CASCADE, related_name="ablocks")
    """The space to which the block is attached."""

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(self, "space", *args, **kwargs)


class BranchBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a branch."""

    class Meta:
        app_label = "lamindb"

    branch: Branch = ForeignKey(Branch, CASCADE, related_name="ablocks")
    """The branch to which the block is attached."""

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(self, "branch", *args, **kwargs)


class ULabelBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a ulabel."""

    class Meta:
        app_label = "lamindb"

    ulabel = ForeignKey("ULabel", CASCADE, related_name="ablocks")
    """The ulabel to which the block is attached."""

    def __init__(self, *args, **kwargs):
        _init_versioned_attached_block(self, "ulabel", *args, **kwargs)
