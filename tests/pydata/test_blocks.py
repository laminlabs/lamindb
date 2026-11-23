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


def test_readme_md_key_is_allowed_and_revises():
    block1 = ln.models.Block(
        key="README.md", content="# v1\n\nhello", kind="readme"
    ).save()
    block2 = ln.models.Block(key="README.md", content="# v2\n\nhello", kind="readme")
    assert block2.stem_uid == block1.stem_uid
    assert block2.uid != block1.uid
    block2.save()
    block1.refresh_from_db()
    assert not block1.is_latest
    block2.delete()
    block1.delete()


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


def test_record_block_readme_same_hash_returns_existing_version():
    """Same readme content on a record returns the latest block."""
    record = ln.Record(name="test-record-blocks").save()
    block1 = ln.models.RecordBlock(record=record, content="1", kind="readme").save()
    block2 = ln.models.RecordBlock(record=record, content="1", kind="readme")
    assert block1 == block2
    assert block1.uid == block2.uid
    changed = ln.models.RecordBlock(record=record, content="2", kind="readme")
    assert changed.stem_uid == block1.stem_uid
    assert changed.uid != block1.uid
    other = ln.Record(name="test-record-blocks-other").save()
    other_block = ln.models.RecordBlock(record=other, content="1", kind="readme")
    assert other_block.uid != block1.uid
    block1.delete()  # BaseSQLRecord has no soft delete; this is permanent
    block2 = ln.models.RecordBlock(record=record, content="1", kind="readme")
    assert block1 != block2  # block2 is a new block (block1 was removed)
    other.delete(permanent=True)
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
    assert "Notes" in result
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

    # same content as the latest readme returns that block
    block_r3.save()  # so next readme for this record gets revises=block_r3
    block_same = ln.models.RecordBlock(record=record, content="v3", kind="readme")
    assert block_same.uid == block_r3.uid
    assert block_same.id == block_r3.id

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

    record.delete(permanent=True)


def test_record_block_filter_respects_default_branch_scope():
    main_branch = ln.Branch.get(name="main")
    ln.setup.switch(main_branch.name)

    main_record = ln.Record(name="record-block-main").save()
    ln.models.RecordBlock(
        record=main_record,
        content="record-block-main-content",
        kind="readme",
        branch=main_branch,
        created_on=main_branch,
    ).save()

    contrib = ln.Branch(name="record_block_scope_branch").save()
    ln.setup.switch(contrib.name)
    contrib_record = ln.Record(name="record-block-contrib").save()
    contrib_block = ln.models.RecordBlock(
        record=contrib_record,
        content="record-block-contrib-content",
        kind="readme",
        branch=contrib,
        created_on=contrib,
    ).save()

    assert (
        ln.models.RecordBlock.filter(content="record-block-contrib-content").count()
        == 1
    )

    ln.setup.switch(main_branch.name)
    assert (
        ln.models.RecordBlock.filter(content="record-block-contrib-content").count()
        == 0
    )

    contrib_block.delete()
    contrib_record.delete(permanent=True)
    main_record.delete(permanent=True)
    contrib.delete(permanent=True)


def test_sqlrecord_notes_returns_latest_readme_content():
    record = ln.Record(name="record-notes-test").save()
    assert record.notes is None

    readme_v1 = ln.models.RecordBlock(
        record=record, content="record notes v1", kind="readme"
    ).save()
    assert record.notes == "record notes v1"

    ln.models.RecordBlock(
        record=record, content="record notes v2", kind="readme"
    ).save()
    readme_v1.refresh_from_db()
    assert not readme_v1.is_latest
    assert record.notes == "record notes v2"

    ln.models.RecordBlock(
        record=record, content="just a comment", kind="comment"
    ).save()
    assert record.notes == "record notes v2"

    record.delete(permanent=True)


def test_branch_notes_returns_latest_readme_content():
    branch = ln.Branch(name="branch-notes-test").save()
    assert branch.notes is None

    readme_v1 = ln.models.BranchBlock(
        branch=branch, content="branch notes v1", kind="readme"
    ).save()
    assert branch.notes == "branch notes v1"

    ln.models.BranchBlock(
        branch=branch, content="branch notes v2", kind="readme"
    ).save()
    readme_v1.refresh_from_db()
    assert not readme_v1.is_latest
    assert branch.notes == "branch notes v2"

    ln.models.BranchBlock(
        branch=branch, content="branch comment", kind="comment"
    ).save()
    assert branch.notes == "branch notes v2"

    branch.emoji = " 🧬 "
    assert branch.emoji == "🧬"
    branch.emoji = None
    assert branch.emoji is None

    branch.delete(permanent=True)


def test_recordblock_revision_ancestry_across_branches():
    from django.db.models.deletion import RestrictedError

    record = ln.Record(name="note-ancestry").save()
    branch = ln.Branch(name="note-ancestry-branch").save()
    sibling = ln.Branch(name="note-ancestry-sibling").save()
    try:
        initial = ln.models.RecordBlock(
            record=record, content="main", kind="readme"
        ).save()
        other = ln.models.RecordBlock(
            record=record, content="sibling", kind="readme", branch=sibling
        ).save()
        first = ln.models.RecordBlock(
            record=record, content="branch first", kind="readme", branch=branch
        ).save()
        assert first.previous_version_id == initial.id
        assert other.previous_version_id == initial.id
        assert ln.models.RecordBlock.objects.get(id=initial.id).is_latest
        second = ln.models.RecordBlock(
            record=record, content="branch second", kind="readme", branch=branch
        ).save()
        assert second.previous_version_id == first.id
        first.refresh_from_db()
        assert first.previous_version_id == initial.id
        assert not first.is_latest

        # Merge moves rows and reconciles heads, without changing ancestry/origin.
        ln.setup.merge(branch, target="main")
        second.refresh_from_db()
        assert second.branch_id == 1
        assert second.created_on_id == branch.id
        assert second.previous_version_id == first.id
        later = ln.models.RecordBlock(
            record=record, content="later main", kind="readme"
        ).save()
        assert later.previous_version_id == second.id
        second.refresh_from_db()
        assert not second.is_latest
        assert second.created_on_id == branch.id

        with pytest.raises(RestrictedError):
            initial.delete()
    finally:
        # RESTRICT permits removal of the entire parent-owned history.
        record.delete(permanent=True)
        branch.delete(permanent=True)
        sibling.delete(permanent=True)


def test_recordblock_records_refreshed_predecessor():
    record = ln.Record(name="note-refreshed-predecessor").save()
    try:
        initial = ln.models.RecordBlock(
            record=record, content="initial", kind="readme"
        ).save()
        pending = ln.models.RecordBlock(record=record, content="pending", kind="readme")
        concurrent = ln.models.RecordBlock(
            record=record, content="concurrent", kind="readme"
        ).save()
        pending.save()
        assert pending.previous_version_id == concurrent.id
        assert concurrent.previous_version_id == initial.id
        assert pending.uid != concurrent.uid
    finally:
        record.delete(permanent=True)


def test_recordblock_run_is_creation_only(monkeypatch):
    record = ln.Record(name="note-producing-run").save()
    empty_record = ln.Record(name="empty-note-producing-run").save()
    transform = ln.Transform(key="note-producing-run.py").save()
    run = ln.Run(transform).save()
    other_run = ln.Run(transform).save()
    try:
        legacy = ln.models.RecordBlock(
            record=record, content="legacy", kind="readme", run=None
        ).save()
        blank = ln.models.RecordBlock(
            record=empty_record, kind="readme", run=None
        ).save()
        monkeypatch.setattr(ln.context, "_run", run)
        reused = ln.models.RecordBlock(
            record=record, content="legacy", kind="readme"
        ).save()
        assert reused.id == legacy.id
        assert reused.run_id is None
        reused_blank = ln.models.RecordBlock(record=empty_record, kind="readme").save()
        assert reused_blank.id == blank.id
        assert reused_blank.run_id is None

        output = ln.models.RecordBlock(
            record=record, content="new", kind="readme"
        ).save()
        assert output.run_id == run.id
        assert output.previous_version_id == legacy.id
        filled = ln.models.RecordBlock(
            record=empty_record, content="filled", kind="readme"
        ).save()
        assert filled.id != blank.id
        assert filled.previous_version_id == blank.id
        assert filled.run_id == run.id
        assert set(run.output_record_blocks.values_list("id", flat=True)) == {
            output.id,
            filled.id,
        }

        monkeypatch.setattr(ln.context, "_run", other_run)
        reused_output = ln.models.RecordBlock(
            record=record, content="new", kind="readme"
        ).save()
        assert reused_output.id == output.id
        assert reused_output.run_id == run.id
        loaded = ln.models.RecordBlock.objects.get(id=output.id)
        loaded.save()
        assert loaded.run_id == run.id
        untracked = ln.models.RecordBlock(
            record=record, content="explicit untracked", kind="readme", run=None
        ).save()
        assert untracked.run_id is None
        explicit = ln.models.RecordBlock(
            record=record, content="explicit run", kind="readme", run=run
        ).save()
        assert explicit.run_id == run.id
        assert not other_run.output_record_blocks.exists()
    finally:
        monkeypatch.setattr(ln.context, "_run", None)
        record.delete(permanent=True)
        empty_record.delete(permanent=True)
        run.delete(permanent=True)
        other_run.delete(permanent=True)
        transform.delete(permanent=True)


def test_recordblock_uid_collision_preserves_existing_version():
    from django.db import IntegrityError

    record = ln.Record(name="note-collision").save()
    branch = ln.Branch(name="note-collision-branch").save()
    sibling = ln.Branch(name="note-collision-sibling").save()
    transform = ln.Transform(key="note-collision.py").save()
    run = ln.Run(transform).save()
    other_run = ln.Run(transform).save()
    try:
        initial = ln.models.RecordBlock(
            record=record, content="main", kind="readme", run=None
        ).save()
        first = ln.models.RecordBlock(
            record=record, content="branch", kind="readme", branch=branch, run=None
        ).save()
        pending = ln.models.RecordBlock(
            record=record, content="pending", kind="readme", branch=branch, run=run
        )
        other = ln.models.RecordBlock(
            record=record,
            content="sibling",
            kind="readme",
            branch=sibling,
            run=other_run,
        ).save()
        # The destination head is unchanged, but another branch used the pending UID.
        with pytest.raises(IntegrityError):
            pending.save()
        other.refresh_from_db()
        assert other.content == "sibling"
        assert other.branch_id == sibling.id
        assert other.created_on_id == sibling.id
        assert other.run_id == other_run.id
        assert other.previous_version_id == initial.id
        first.refresh_from_db()
        assert first.is_latest
        assert not run.output_record_blocks.exists()
        assert record.ablocks.count() == 3
    finally:
        record.delete(permanent=True)
        branch.delete(permanent=True)
        sibling.delete(permanent=True)
        run.delete(permanent=True)
        other_run.delete(permanent=True)
        transform.delete(permanent=True)


def test_recordblock_rejects_unrelated_predecessor():
    first = ln.Record(name="note-owner-first").save()
    second = ln.Record(name="note-owner-second").save()
    try:
        original = ln.models.RecordBlock(
            record=first, content="initial", kind="readme"
        ).save()
        with pytest.raises(ValueError, match="same record and kind"):
            ln.models.RecordBlock(
                record=second, content="changed", kind="readme", revises=original
            )
        comment = ln.models.RecordBlock(
            record=first, content="comment", kind="comment"
        ).save()
        with pytest.raises(ValueError, match="same record and kind"):
            ln.models.RecordBlock(
                record=first, content="changed", kind="readme", revises=comment
            )
        with pytest.raises(ValueError, match="version family"):
            ln.models.RecordBlock(
                record=first,
                content="changed",
                kind="readme",
                revises=original,
                uid="otherfamily000000001",
            )
    finally:
        first.delete(permanent=True)
        second.delete(permanent=True)
