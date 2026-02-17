from __future__ import annotations

from typing import TYPE_CHECKING, Any, overload

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

    from ..base.types import BlockKind
    from .query_manager import RelatedManager

_VERSIONED_ATTACHED_KINDS = ("readme",)  # only readme is versioned; comment is not
_VALID_BLOCK_KINDS: tuple[str, ...] = ("readme", "comment")


def _init_versioned_attached_block(
    self: BaseBlock,
    fk_field_name: str,
    fk_value: Any,
    content: str | None,
    kind: str,
    version_tag: str | None,
    revises: IsVersioned | None,
    uid: str | None,
    skip_hash_lookup: bool,
    using_key: str | None,
    **extra_kwargs: Any,
) -> None:
    cls = type(self)
    if kind == "comment" and revises is not None:
        raise ValueError(
            "revises is not allowed for kind='comment'; comments are not versioned"
        )
    if kind == "readme" and fk_value is not None:
        if revises is None:
            candidate_for_revises = (
                cls.objects.using(using_key)
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
        return None
    if kind == "comment":
        new_uid, revises = create_uid(
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
    if revises is not None or uid is not None:
        new_uid, revises = create_uid(
            revises=revises,
            version_tag=version_tag,
            n_full_id=cls._len_full_uid,
        )
        if uid is None:
            uid = new_uid
    else:
        new_uid, revises = create_uid(
            revises=None,
            version_tag=version_tag,
            n_full_id=cls._len_full_uid,
        )
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
    """A markdown block for anything: issues, standalone markdown pages, comments, etc."""

    class Meta:
        app_label = "lamindb"
        unique_together = ("key", "hash")

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
        kind: BlockKind = ...,
        version: str | None = None,
        revises: Block | None = None,
        anchor: Block | None = None,
        skip_hash_lookup: bool = False,
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
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        branch = kwargs.pop("branch", None)
        branch_id = kwargs.pop("branch_id", 1)
        space = kwargs.pop("space", None)
        space_id = kwargs.pop("space_id", 1)
        if kwargs:
            raise ValueError(
                "Only key, content, kind, version, revises, anchor, "
                f"skip_hash_lookup can be passed, but you passed: {kwargs}"
            )
        if kind is None:
            raise ValueError("kind is required for Block; use 'readme' or 'comment'")
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        if revises is None:
            if uid is not None:
                revises = (
                    Block.objects.using(using_key)
                    .filter(
                        uid__startswith=uid[:-4],
                        is_latest=True,
                    )
                    .order_by("-created_at")
                    .first()
                )
            elif key is not None:
                candidate_for_revises = (
                    Block.objects.using(using_key)
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
        if content is not None and not skip_hash_lookup:
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
                        f"passed key {key}, keeping original key; update "
                        "manually if needed or pass skip_hash_lookup if you "
                        "want to duplicate the block"
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
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a RecordBlock"
            )
        record = kwargs.pop("record", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        if kwargs:
            raise ValueError(
                "Only record, content, kind, version, revises, skip_hash_lookup "
                f"can be passed, but you passed: {kwargs}"
            )
        if record is None:
            raise ValueError("record is required for RecordBlock")
        if kind is None:
            raise ValueError(
                "kind is required for RecordBlock; use 'readme' or 'comment'"
            )
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "record",
            record,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
        )


class ArtifactBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to an artifact."""

    class Meta:
        app_label = "lamindb"

    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="ablocks")
    """The artifact to which the block is attached."""

    def __init__(self, *args, **kwargs):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct an ArtifactBlock"
            )
        artifact = kwargs.pop("artifact", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        if kwargs:
            raise ValueError(
                "Only artifact, content, kind, version, revises, skip_hash_lookup "
                f"can be passed, but you passed: {kwargs}"
            )
        if artifact is None:
            raise ValueError("artifact is required for ArtifactBlock")
        if kind is None:
            raise ValueError(
                "kind is required for ArtifactBlock; use 'readme' or 'comment'"
            )
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "artifact",
            artifact,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
        )


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
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a TransformBlock"
            )
        transform = kwargs.pop("transform", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        line_number = kwargs.pop("line_number", None)
        if kwargs:
            raise ValueError(
                "Only transform, content, kind, version, revises, skip_hash_lookup, "
                f"line_number can be passed, but you passed: {kwargs}"
            )
        if transform is None:
            raise ValueError("transform is required for TransformBlock")
        if kind is None:
            raise ValueError(
                "kind is required for TransformBlock; use 'readme' or 'comment'"
            )
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "transform",
            transform,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
            line_number=line_number,
        )


class RunBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a run."""

    class Meta:
        app_label = "lamindb"

    run: Run = ForeignKey(Run, CASCADE, related_name="ablocks")
    """The run to which the block is attached."""

    def __init__(self, *args, **kwargs):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a RunBlock"
            )
        run = kwargs.pop("run", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        if kwargs:
            raise ValueError(
                "Only run, content, kind, version, revises, skip_hash_lookup "
                f"can be passed, but you passed: {kwargs}"
            )
        if run is None:
            raise ValueError("run is required for RunBlock")
        if kind is None:
            raise ValueError("kind is required for RunBlock; use 'readme' or 'comment'")
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "run",
            run,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
        )


class CollectionBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a collection."""

    class Meta:
        app_label = "lamindb"

    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="ablocks", null=True
    )
    """The collection to which the block is attached."""

    def __init__(self, *args, **kwargs):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a CollectionBlock"
            )
        collection = kwargs.pop("collection", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        if kwargs:
            raise ValueError(
                "Only collection, content, kind, version, revises, skip_hash_lookup "
                f"can be passed, but you passed: {kwargs}"
            )
        if collection is None:
            raise ValueError("collection is required for CollectionBlock")
        if kind is None:
            raise ValueError(
                "kind is required for CollectionBlock; use 'readme' or 'comment'"
            )
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "collection",
            collection,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
        )


class SchemaBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a schema."""

    class Meta:
        app_label = "lamindb"

    schema: Schema = ForeignKey(Schema, CASCADE, related_name="ablocks")
    """The schema to which the block is attached."""

    def __init__(self, *args, **kwargs):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a SchemaBlock"
            )
        schema = kwargs.pop("schema", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        if kwargs:
            raise ValueError(
                "Only schema, content, kind, version, revises, skip_hash_lookup "
                f"can be passed, but you passed: {kwargs}"
            )
        if schema is None:
            raise ValueError("schema is required for SchemaBlock")
        if kind is None:
            raise ValueError(
                "kind is required for SchemaBlock; use 'readme' or 'comment'"
            )
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "schema",
            schema,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
        )


class FeatureBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a feature."""

    class Meta:
        app_label = "lamindb"

    feature: Feature = ForeignKey(Feature, CASCADE, related_name="ablocks")
    """The feature to which the block is attached."""

    def __init__(self, *args, **kwargs):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a FeatureBlock"
            )
        feature = kwargs.pop("feature", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        if kwargs:
            raise ValueError(
                "Only feature, content, kind, version, revises, skip_hash_lookup "
                f"can be passed, but you passed: {kwargs}"
            )
        if feature is None:
            raise ValueError("feature is required for FeatureBlock")
        if kind is None:
            raise ValueError(
                "kind is required for FeatureBlock; use 'readme' or 'comment'"
            )
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "feature",
            feature,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
        )


class ProjectBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a project."""

    class Meta:
        app_label = "lamindb"

    project: Project = ForeignKey(Project, CASCADE, related_name="ablocks")
    """The project to which the block is attached."""

    def __init__(self, *args, **kwargs):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a ProjectBlock"
            )
        project = kwargs.pop("project", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        if kwargs:
            raise ValueError(
                "Only project, content, kind, version, revises, skip_hash_lookup "
                f"can be passed, but you passed: {kwargs}"
            )
        if project is None:
            raise ValueError("project is required for ProjectBlock")
        if kind is None:
            raise ValueError(
                "kind is required for ProjectBlock; use 'readme' or 'comment'"
            )
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "project",
            project,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
        )


class SpaceBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a space."""

    class Meta:
        app_label = "lamindb"

    space: Space = ForeignKey(Space, CASCADE, related_name="ablocks")
    """The space to which the block is attached."""

    def __init__(self, *args, **kwargs):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a SpaceBlock"
            )
        space = kwargs.pop("space", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        if kwargs:
            raise ValueError(
                "Only space, content, kind, version, revises, skip_hash_lookup "
                f"can be passed, but you passed: {kwargs}"
            )
        if space is None:
            raise ValueError("space is required for SpaceBlock")
        if kind is None:
            raise ValueError(
                "kind is required for SpaceBlock; use 'readme' or 'comment'"
            )
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "space",
            space,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
        )


class BranchBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a branch."""

    class Meta:
        app_label = "lamindb"

    branch: Branch = ForeignKey(Branch, CASCADE, related_name="ablocks")
    """The branch to which the block is attached."""

    def __init__(self, *args, **kwargs):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a BranchBlock"
            )
        branch = kwargs.pop("branch", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        if kwargs:
            raise ValueError(
                "Only branch, content, kind, version, revises, skip_hash_lookup "
                f"can be passed, but you passed: {kwargs}"
            )
        if branch is None:
            raise ValueError("branch is required for BranchBlock")
        if kind is None:
            raise ValueError(
                "kind is required for BranchBlock; use 'readme' or 'comment'"
            )
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "branch",
            branch,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
        )


class ULabelBlock(BaseBlock, BaseSQLRecord):
    """An unstructured notes block that can be attached to a ulabel."""

    class Meta:
        app_label = "lamindb"

    ulabel = ForeignKey("ULabel", CASCADE, related_name="ablocks")
    """The ulabel to which the block is attached."""

    def __init__(self, *args, **kwargs):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if args:
            raise ValueError(
                "Please only use keyword arguments to construct a ULabelBlock"
            )
        ulabel = kwargs.pop("ulabel", None)
        content = kwargs.pop("content", None)
        kind = kwargs.pop("kind", None)
        version_tag = kwargs.pop("version_tag", kwargs.pop("version", None))
        revises = kwargs.pop("revises", None)
        skip_hash_lookup = kwargs.pop("skip_hash_lookup", False)
        using_key = kwargs.pop("using_key", None)
        uid = kwargs.pop("uid", None) if "uid" in kwargs else None
        if kwargs:
            raise ValueError(
                "Only ulabel, content, kind, version, revises, skip_hash_lookup "
                f"can be passed, but you passed: {kwargs}"
            )
        if ulabel is None:
            raise ValueError("ulabel is required for ULabelBlock")
        if kind is None:
            raise ValueError(
                "kind is required for ULabelBlock; use 'readme' or 'comment'"
            )
        if kind not in _VALID_BLOCK_KINDS:
            raise ValueError(f"kind must be 'readme' or 'comment', got {kind!r}")
        _init_versioned_attached_block(
            self,
            "ulabel",
            ulabel,
            content,
            kind,
            version_tag,
            revises,
            uid,
            skip_hash_lookup,
            using_key,
        )
