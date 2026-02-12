"""Tests for ln.setup.merge."""

import lamindb as ln
import pytest


def test_merge_branch_into_main():
    """Merge a branch into main: create branch, add ULabel, switch to main, merge."""
    branch = ln.Branch(name="test_merge_branch").save()
    ln.setup.switch(branch.name)
    assert ln.setup.settings.branch == branch
    assert ln.setup.settings.branch.name == "test_merge_branch"

    ulabel = ln.ULabel(name="test_merge_record").save()
    assert ulabel.branch == branch

    ln.setup.switch("main")
    assert ln.setup.settings.branch.name == "main"
    assert ln.ULabel.filter(name="test_merge_record").count() == 0

    ln.setup.merge("test_merge_branch")
    assert ln.ULabel.filter(name="test_merge_record").count() == 1
    ulabel = ln.ULabel.get(name="test_merge_record")
    assert ulabel.branch.name == "main"

    # Clean up
    ulabel.delete(permanent=True)
    branch.delete(permanent=True)
    ln.setup.switch("main")


def test_merge_nonexistent_branch_raises():
    """Merge a non-existent branch raises ObjectDoesNotExist."""
    with pytest.raises(ln.errors.ObjectDoesNotExist) as exc_info:
        ln.setup.merge("nonexistent_branch_xyz")
    assert "not found" in str(exc_info.value).lower()
