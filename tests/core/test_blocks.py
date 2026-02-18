import lamindb as ln
import pytest


def test_block_recovery_based_on_hash():
    block1 = ln.models.Block(key="__lamindb_block__", content="1", kind="readme").save()
    block2 = ln.models.Block(key="__lamindb_block__", content="1", kind="readme")
    assert block1 == block2
    block1.delete()
    block2 = ln.models.Block(key="__lamindb_block__", content="1", kind="readme")
    assert block1 != block2
    block1.delete(permanent=True)


def test_block_recovery_based_on_key():
    block1 = ln.models.Block(key="__lamindb_block__", kind="readme").save()
    block2 = ln.models.Block(key="__lamindb_block__", kind="readme")
    assert block1 == block2
    block1.delete()
    block2 = ln.models.Block(key="__lamindb_block__", kind="readme")
    assert block1 != block2
    block1.delete(permanent=True)


def test_revise_blocks():
    # attempt to create a block with an invalid version
    with pytest.raises(ValueError) as error:
        ln.models.Block(key="__lamindb_block__", version=0, kind="readme")
    assert "version" in error.exconly() or "version_tag" in error.exconly()

    # create a versioned block
    block = ln.models.Block(key="__lamindb_block__", version="1", kind="readme")
    assert block.version_tag == "1"
    assert block.version == "1"
    assert len(block.uid) == ln.models.Block._len_full_uid == 20
    assert len(block.stem_uid) == ln.models.Block._len_stem_uid == 16

    block.save()

    # try to reload the same block with the same uid
    block_reload = ln.models.Block(
        uid=block.uid, key="__lamindb_artifact__", kind="readme"
    )
    assert block_reload.id == block.id
    assert block_reload.key == "__lamindb_block__"  # unchanged, prints logging

    # create new block from old block
    block_r2 = ln.models.Block(content="v2", revises=block, kind="readme")
    assert block_r2.uid != block.uid
    assert block_r2.uid.endswith("0001")
    block_r2 = ln.models.Block(content="v2", revises=block, kind="readme")
    assert block_r2.uid != block.uid
    assert block_r2.uid.endswith("0001")
    assert block_r2.stem_uid == block.stem_uid
    assert block_r2.version_tag is None
    assert block_r2.version == block_r2.uid[-4:]
    assert block_r2.is_latest
    assert block.is_latest
    block_r2.save()
    assert not block.is_latest

    # create new block from newly versioned block
    block_r3 = ln.models.Block(
        content="v3", revises=block_r2, version="2", kind="readme"
    )
    assert block_r3.stem_uid == block.stem_uid
    assert block_r3.version_tag == "2"
    assert block_r3.version == "2"

    # revise by matching on key
    key = "__lamindb_artifact__"
    block_r2.key = key
    block_r2.save()
    assert block_r2.is_latest
    block_r3 = ln.models.Block(content="v3", key=key, version="2", kind="readme")
    assert block_r3.uid[:-4] == block_r2.uid[:-4]
    assert block_r3.uid != block_r2.uid  # new version after block_r2
    block_r2.content = "something else"
    block_r2.save()
    block_r3 = ln.models.Block(content="v3", key=key, version="2", kind="readme")
    assert block_r3.uid[:-4] == block_r2.uid[:-4]
    assert block_r3.uid != block_r2.uid  # yet another new version
    assert block_r3.stem_uid == block_r2.stem_uid
    assert block_r3.key == key
    assert block_r3.version_tag == "2"
    assert block_r3.version == "2"
    assert block_r3.is_latest
    assert block_r2.is_latest
    assert block_r3._revises is not None
    block_r3.save()
    block_r2 = ln.models.Block.get(block_r2.uid)
    assert not block_r2.is_latest

    # wrong block type
    with pytest.raises(TypeError) as error:
        ln.models.Block(
            key="__lamindb_block__", revises=ln.Record(name="x"), kind="readme"
        )
    assert error.exconly().startswith("TypeError: `revises` has to be of type `Block`")

    # wrong kwargs
    with pytest.raises(ValueError) as error:
        ln.models.Block(key="__lamindb_block__", x=1, kind="readme")
    assert "can be passed" in error.exconly() and "x" in error.exconly()

    # kind required (Block only supports kind="readme")
    with pytest.raises(ValueError) as error:
        ln.models.Block(key="__lamindb_block__", content="y")
    assert "kind" in error.exconly() and "readme" in error.exconly()

    # invalid kind (Block only supports readme)
    with pytest.raises(ValueError) as error:
        ln.models.Block(key="__lamindb_block__", content="y", kind="comment")
    assert "readme" in error.exconly() or "Only kind" in error.exconly()

    # cleanup
    block_r2.delete()
    block.delete()

    # unversioned block
    block = ln.models.Block(key="__lamindb_block__", kind="readme")
    assert block.version_tag is None
    assert block.version == block.uid[-4:]
    block.save()

    # create new block from old block
    new_block = ln.models.Block(content="new", revises=block, kind="readme")
    assert block.version_tag is None
    assert block.version == block.uid[-4:]
    assert new_block.stem_uid == block.stem_uid
    assert new_block.uid.endswith("0001")
    assert new_block.version_tag is None
    assert new_block.version == new_block.uid[-4:]

    block.delete(permanent=True)


def test_record_block_readme_always_new_version():
    """Readme always creates a new version (no content-hash dedup)."""
    record = ln.Record(name="test-record-blocks").save()
    block1 = ln.models.RecordBlock(record=record, content="1", kind="readme").save()
    block2 = ln.models.RecordBlock(record=record, content="1", kind="readme")
    assert block1.stem_uid == block2.stem_uid
    assert block1.uid != block2.uid  # new version each time
    block1.delete()  # BaseSQLRecord has no soft delete; this is permanent
    block2 = ln.models.RecordBlock(record=record, content="1", kind="readme")
    assert block1 != block2  # block2 is a new block (block1 was removed)
    record.delete(permanent=True)


def test_record_block_comment_always_new_block():
    """Comment always creates a new block (no versioning; revises not allowed)."""
    record = ln.Record(name="test-record-blocks-comment").save()
    # Add readme and comments to test full describe
    ln.models.RecordBlock(
        record=record, content="# Overview\n\nTest readme.", kind="readme"
    ).save()
    # Comments never version: each creation is a new comment (new uid).
    comment1 = ln.models.RecordBlock(
        record=record, content="same text", kind="comment"
    ).save()
    comment2 = ln.models.RecordBlock(record=record, content="same text", kind="comment")
    assert comment1.stem_uid != comment2.stem_uid  # always new comment, no dedup
    # revises is not allowed for kind='comment'
    with pytest.raises(ValueError) as error:
        ln.models.RecordBlock(
            record=record, content="a comment", kind="comment", revises=comment1
        )
    assert "revises is not allowed for kind='comment'" in error.exconly()

    # Test full describe call with include="comments"
    result = record.describe(return_str=True, include="comments")
    assert "README" in result
    assert "comment by" in result
    assert "same text" in result

    comment1.delete()
    record.delete(permanent=True)


def test_record_block_recovery_based_on_record_and_kind():
    record = ln.Record(name="test-record-blocks-key").save()
    block1 = ln.models.RecordBlock(record=record, kind="readme").save()
    block2 = ln.models.RecordBlock(record=record, kind="readme")
    assert block1 == block2
    block1.delete()  # BaseSQLRecord has no soft delete; this is permanent
    block2 = ln.models.RecordBlock(record=record, kind="readme")
    assert block1 != block2  # block2 is a new block (block1 was removed)
    record.delete(permanent=True)


def test_revise_record_blocks():
    record = ln.Record(name="test-record-revise").save()

    # create a versioned record block
    block = ln.models.RecordBlock(
        record=record, content="v1", kind="readme", version="1"
    )
    assert block.version_tag == "1"
    assert block.version == "1"
    assert len(block.uid) == ln.models.RecordBlock._len_full_uid == 20
    assert len(block.stem_uid) == ln.models.RecordBlock._len_stem_uid == 16
    block.save()

    # reload same block by uid
    block_reload = ln.models.RecordBlock(record=record, uid=block.uid, kind="readme")
    assert block_reload.id == block.id

    # create new block from old block
    block_r2 = ln.models.RecordBlock(
        record=record, content="v2", kind="readme", revises=block
    )
    assert block_r2.uid != block.uid
    assert block_r2.uid.endswith("0001")
    assert block_r2.stem_uid == block.stem_uid
    assert block_r2.is_latest
    assert block.is_latest
    block_r2.save()
    assert not block.is_latest

    # create new block from newly versioned block
    block_r3 = ln.models.RecordBlock(
        record=record, content="v3", kind="readme", revises=block_r2, version="2"
    )
    assert block_r3.stem_uid == block.stem_uid
    assert block_r3.version_tag == "2"
    assert block_r3.version == "2"

    # readme always creates a new version (no hash-based dedup)
    block_r3.save()  # so next readme for this record gets revises=block_r3
    block_same = ln.models.RecordBlock(record=record, content="v3", kind="readme")
    assert block_same.stem_uid == block_r3.stem_uid
    assert block_same.uid != block_r3.uid  # new version (0003)

    # comment does not accept revises
    with pytest.raises(ValueError) as error:
        ln.models.RecordBlock(
            record=record, content="a comment", kind="comment", revises=block
        )
    assert "revises is not allowed for kind='comment'" in error.exconly()

    # wrong kwargs
    with pytest.raises(ValueError) as error:
        ln.models.RecordBlock(record=record, x=1)
    assert "can be passed" in error.exconly()

    # record required
    with pytest.raises(ValueError) as error:
        ln.models.RecordBlock(content="x", kind="readme")
    assert "record is required" in error.exconly()

    block_r2.delete()
    block.delete()
    record.delete(permanent=True)
