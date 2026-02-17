import lamindb as ln
import pytest


def test_block_recovery_based_on_hash():
    block1 = ln.models.Block(key="my-block", content="1").save()
    block2 = ln.models.Block(key="my-block", content="1")
    assert block1 == block2
    block1.delete()
    block2 = ln.models.Block(key="my-block", content="1")
    assert block1 != block2
    block1.delete(permanent=True)


def test_block_recovery_based_on_key():
    block1 = ln.models.Block(key="my-block").save()
    block2 = ln.models.Block(key="my-block")
    assert block1 == block2
    block1.delete()
    block2 = ln.models.Block(key="my-block")
    assert block1 != block2
    block1.delete(permanent=True)


def test_revise_blocks():
    # attempt to create a block with an invalid version
    with pytest.raises(ValueError) as error:
        ln.models.Block(key="My block", version=0)
    assert "version" in error.exconly() or "version_tag" in error.exconly()

    # create a versioned block
    block = ln.models.Block(key="My block", version="1")
    assert block.version_tag == "1"
    assert block.version == "1"
    assert len(block.uid) == ln.models.Block._len_full_uid == 20
    assert len(block.stem_uid) == ln.models.Block._len_stem_uid == 16

    block.save()

    # try to reload the same block with the same uid
    block_reload = ln.models.Block(uid=block.uid, key="My block updated name")
    assert block_reload.id == block.id
    assert block_reload.key == "My block"  # unchanged, prints logging

    # create new block from old block
    block_r2 = ln.models.Block(content="v2", revises=block)
    assert block_r2.uid != block.uid
    assert block_r2.uid.endswith("0001")
    block_r2 = ln.models.Block(content="v2", revises=block)
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
    block_r3 = ln.models.Block(content="v3", revises=block_r2, version="2")
    assert block_r3.stem_uid == block.stem_uid
    assert block_r3.version_tag == "2"
    assert block_r3.version == "2"

    # revise by matching on key
    key = "my-readme.md"
    block_r2.key = key
    block_r2.save()
    assert block_r2.is_latest
    block_r3 = ln.models.Block(content="v3", key=key, version="2")
    assert block_r3.uid[:-4] == block_r2.uid[:-4]
    assert block_r3.uid != block_r2.uid  # new version after block_r2
    block_r2.content = "something else"
    block_r2.save()
    block_r3 = ln.models.Block(content="v3", key=key, version="2")
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
        ln.models.Block(revises=ln.Record(name="x"))
    assert error.exconly().startswith("TypeError: `revises` has to be of type `Block`")

    # wrong kwargs
    with pytest.raises(ValueError) as error:
        ln.models.Block(x=1)
    assert "can be passed" in error.exconly() and "x" in error.exconly()

    # cleanup
    block_r2.delete()
    block.delete()

    # unversioned block
    block = ln.models.Block(key="My block")
    assert block.version_tag is None
    assert block.version == block.uid[-4:]
    block.save()

    # create new block from old block
    new_block = ln.models.Block(content="new", revises=block)
    assert block.version_tag is None
    assert block.version == block.uid[-4:]
    assert new_block.stem_uid == block.stem_uid
    assert new_block.uid.endswith("0001")
    assert new_block.version_tag is None
    assert new_block.version == new_block.uid[-4:]

    block.delete(permanent=True)


def test_record_block_recovery_based_on_hash():
    record = ln.Record(name="test-record-blocks").save()
    block1 = ln.models.RecordBlock(record=record, content="1", kind="mdpage").save()
    block2 = ln.models.RecordBlock(record=record, content="1", kind="mdpage")
    assert block1 == block2
    block1.delete()  # BaseSQLRecord has no soft delete; this is permanent
    block2 = ln.models.RecordBlock(record=record, content="1", kind="mdpage")
    assert block1 != block2  # block2 is a new block (block1 was removed)
    record.delete(permanent=True)


def test_record_block_recovery_based_on_record_and_kind():
    record = ln.Record(name="test-record-blocks-key").save()
    block1 = ln.models.RecordBlock(record=record, kind="mdpage").save()
    block2 = ln.models.RecordBlock(record=record, kind="mdpage")
    assert block1 == block2
    block1.delete()  # BaseSQLRecord has no soft delete; this is permanent
    block2 = ln.models.RecordBlock(record=record, kind="mdpage")
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

    # same (record, kind) + same content hash returns existing
    block_same = ln.models.RecordBlock(record=record, content="v3", kind="readme")
    assert block_same.id == block_r3.id

    # wrong kwargs
    with pytest.raises(ValueError) as error:
        ln.models.RecordBlock(record=record, x=1)
    assert "can be passed" in error.exconly()

    # record required
    with pytest.raises(ValueError) as error:
        ln.models.RecordBlock(content="x", kind="mdpage")
    assert "record is required" in error.exconly()

    block_r2.delete()
    block.delete()
    record.delete(permanent=True)
