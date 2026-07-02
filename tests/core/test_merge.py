"""Tests for ln.setup.merge."""

import lamindb as ln
import pytest


def test_merge_branch_into_main():
    """Merge a branch into main: create branch, add ULabel, switch to main, merge."""
    branch = ln.Branch(name="test_merge_branch").save()
    assert branch.status == "standalone"
    ln.setup.switch(branch.name)
    assert ln.setup.settings.branch == branch
    assert ln.setup.settings.branch.name == "test_merge_branch"

    ulabel = ln.ULabel(name="test_merge_record").save()
    assert ulabel.branch == branch
    assert ulabel.created_on == branch  # created_on set to creation branch

    ln.setup.switch("main")
    assert ln.setup.settings.branch.name == "main"
    assert ln.setup.settings.branch.status == "standalone"
    assert ln.ULabel.filter(name="test_merge_record").count() == 0

    ln.setup.merge("test_merge_branch")
    assert ln.ULabel.filter(name="test_merge_record").count() == 1
    ulabel = ln.ULabel.get(name="test_merge_record")
    assert ulabel.branch.name == "main"
    # created_on still points to the branch on which the record was created
    assert ulabel.created_on == branch
    assert ulabel.created_on.name == "test_merge_branch"
    # merged branch has status "merged"
    branch.refresh_from_db()
    assert branch.status == "merged"
    # this is a merge call to check that branch.describe() works because it
    # has a custom describe method
    branch.describe(return_str=True)

    # Clean up
    ulabel.delete(permanent=True)
    branch.delete(permanent=True)
    ln.setup.switch("main")


def test_branch_status_values():
    """Branch status maps codes onto standalone/draft/review/merged/closed."""
    main_branch = ln.Branch.get(name="main")
    assert main_branch.status == "standalone"
    archive_branch = ln.Branch.get(name="archive")
    assert archive_branch.status == "standalone"
    trash_branch = ln.Branch.get(name="trash")
    assert trash_branch.status == "standalone"
    # User-created branch is standalone by default.
    branch = ln.Branch(name="test_status_branch").save()
    assert branch.status == "standalone"
    branch.status = "draft"
    branch.save()
    branch.refresh_from_db()
    assert branch.status == "draft"
    branch.status = "review"
    branch.save()
    branch.refresh_from_db()
    assert branch.status == "review"
    branch.status = "closed"
    branch.save()
    branch.refresh_from_db()
    assert branch.status == "closed"
    branch.delete(permanent=True)


def test_draft_review_and_close_merge_request_status():
    branch = ln.Branch(name="test_mr_draft_review_close").save()
    assert branch.status == "standalone"

    branch.status = "draft"
    branch.save()
    branch.refresh_from_db()
    assert branch.status == "draft"

    branch.status = "review"
    branch.save()
    branch.refresh_from_db()
    assert branch.status == "review"

    branch.status = "closed"
    branch.save()
    branch.refresh_from_db()
    assert branch.status == "closed"

    branch.delete(permanent=True)


def test_merge_nonexistent_branch_raises():
    """Merge a non-existent branch raises ObjectDoesNotExist."""
    with pytest.raises(ln.errors.ObjectDoesNotExist) as exc_info:
        ln.setup.merge("nonexistent_branch_xyz")
    assert "not found" in str(exc_info.value).lower()


def test_merge_without_switch_via_explicit_target():
    main_branch = ln.Branch.get(name="main")
    source_branch = ln.Branch(name="test_merge_explicit_target_source").save()
    ln.setup.switch(source_branch.name)
    ulabel = ln.ULabel(name="test_merge_explicit_target_record").save()
    assert ln.setup.settings.branch.name == source_branch.name
    assert ulabel.branch == source_branch

    ln.setup.merge(source_branch.name, target=main_branch.name)

    ulabel.refresh_from_db()
    assert ln.setup.settings.branch.name == source_branch.name
    assert ulabel.branch.name == main_branch.name
    assert ulabel.created_on == source_branch

    ulabel.delete(permanent=True)
    source_branch.delete(permanent=True)
    ln.setup.switch(main_branch.name)


def test_merge_with_target_branch_object():
    main_branch = ln.Branch.get(name="main")
    source_branch = ln.Branch(name="test_merge_target_object_source").save()
    ln.setup.switch(source_branch.name)
    ulabel = ln.ULabel(name="test_merge_target_object_record").save()
    assert ulabel.branch == source_branch

    ln.setup.merge(source_branch, target=main_branch)

    ulabel.refresh_from_db()
    assert ulabel.branch == main_branch
    assert ulabel.created_on == source_branch

    ulabel.delete(permanent=True)
    source_branch.delete(permanent=True)
    ln.setup.switch(main_branch.name)


def test_merge_explicit_target_not_found_raises():
    source_branch = ln.Branch(name="test_merge_target_not_found_source").save()
    with pytest.raises(ln.errors.ObjectDoesNotExist) as exc_info:
        ln.setup.merge(source_branch.name, target="nonexistent_target_branch_xyz")
    assert "not found" in str(exc_info.value).lower()
    source_branch.delete(permanent=True)


def test_merge_reconciles_is_latest_for_versioned_records():
    main_branch = ln.Branch.get(name="main")
    ln.setup.switch(main_branch.name)

    transform_v1 = ln.Transform(
        key="test-merge-is-latest",
        source_code="main-v1",
        kind="pipeline",
    ).save()
    branch = ln.Branch(name="test_merge_latest_branch").save()
    ln.setup.switch(branch.name)
    transform_v2 = ln.Transform(
        key="test-merge-is-latest",
        revises=transform_v1,
        source_code="feature-v2",
        kind="pipeline",
    ).save()
    transform_v1.refresh_from_db()
    assert transform_v1.is_latest
    assert transform_v2.is_latest

    ln.setup.switch(main_branch.name)
    ln.setup.merge(branch.name)

    family = ln.Transform.objects.filter(
        uid__startswith=transform_v1.uid[:-4], branch_id=1
    )
    assert family.filter(is_latest=True).count() == 1
    assert family.get(is_latest=True).uid == transform_v2.uid

    for record in family:
        record.delete(permanent=True)
    branch.delete(permanent=True)


def test_merge_updates_recordblock_branch():
    main_branch = ln.Branch.get(name="main")
    ln.setup.switch(main_branch.name)

    source_branch = ln.Branch(name="test_merge_recordblock_branch").save()
    ln.setup.switch(source_branch.name)
    record = ln.Record(name="recordblock-merge-record").save()
    block = ln.models.RecordBlock(
        record=record,
        content="recordblock merge content",
        kind="readme",
        branch=source_branch,
        created_on=source_branch,
    ).save()
    assert block.branch == source_branch
    assert block.created_on == source_branch

    ln.setup.switch(main_branch.name)
    ln.setup.merge(source_branch.name)

    block.refresh_from_db()
    assert block.branch.name == "main"
    assert block.created_on == source_branch

    record.delete(permanent=True)
    source_branch.delete(permanent=True)
