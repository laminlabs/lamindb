"""Tests for ln.setup.switch."""

from pathlib import Path

import lamindb as ln
import pytest

from worktree_test_utils import run_worktree_case


def test_switch_create_existing_branch_raises():
    """Switch with create=True and existing branch raises BranchAlreadyExists with hint."""
    with pytest.raises(ln.errors.BranchAlreadyExists) as exc_info:
        ln.setup.switch("main", create=True)
    msg = str(exc_info.value)
    assert "already exists" in msg
    assert "-c/--create" in msg or "Omit" in msg


def test_switch_space_does_not_depend_on_worktree_bootstrap_state():
    ln.setup.switch("all", space=True)
    assert ln.setup.settings.space.name == "all"


def test_switch_create_worktree_from_dev_dir_root(tmp_path: Path):
    run_worktree_case("switch_creates_worktree", tmp_path)
