"""Tests for ln.setup.merge."""

import lamindb as ln
import pytest


def test_merge_branch_into_main():
    """Merge a branch into main: create branch, add ULabel, switch to main, merge."""
    branch = ln.Branch(name="test_merge_branch").save()
    assert branch.status == "open"
    ln.setup.switch(branch.name)
    assert ln.setup.settings.branch == branch
    assert ln.setup.settings.branch.name == "test_merge_branch"

    ulabel = ln.ULabel(name="test_merge_record").save()
    assert ulabel.branch == branch
    assert ulabel.created_on == branch  # created_on set to creation branch

    ln.setup.switch("main")
    assert ln.setup.settings.branch.name == "main"
    assert ln.setup.settings.branch.status == "builtin"
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
    """Builtin branches have status 'builtin', new branches 'open'."""
    main_branch = ln.Branch.get(name="main")
    assert main_branch.status == "builtin"
    archive_branch = ln.Branch.get(name="archive")
    assert archive_branch.status == "builtin"
    trash_branch = ln.Branch.get(name="trash")
    assert trash_branch.status == "builtin"
    # User-created branch is "open" until merged
    branch = ln.Branch(name="test_status_branch").save()
    assert branch.status == "open"
    branch.delete(permanent=True)


def test_merge_nonexistent_branch_raises():
    """Merge a non-existent branch raises ObjectDoesNotExist."""
    with pytest.raises(ln.errors.ObjectDoesNotExist) as exc_info:
        ln.setup.merge("nonexistent_branch_xyz")
    assert "not found" in str(exc_info.value).lower()


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
