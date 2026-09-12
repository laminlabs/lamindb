"""Tests for ln.setup.switch."""

import os
import shutil
import time
from pathlib import Path

import lamindb as ln
import lamindb_setup as ln_setup
import pytest
from lamindb_setup.core._settings_store import local_current_branch_file


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
    previous_dev_dir = ln_setup.settings.dev_dir
    previous_worktree = ln_setup.settings.worktree
    previous_cwd = Path.cwd()
    worktree_parent = tmp_path / "worktrees"
    worktree_parent.mkdir(parents=True, exist_ok=True)
    branch_name = f"wt-{time.time_ns()}"
    child = worktree_parent / branch_name
    try:
        ln_setup.settings.dev_dir = worktree_parent
        ln_setup.settings.worktree = True
        os.chdir(worktree_parent)
        ln.setup.switch(branch_name, create=True)
        assert child.exists()
        assert (
            local_current_branch_file(child).read_text().split("\n")[1] == branch_name
        )
    finally:
        os.chdir(previous_cwd)
        ln.Branch.filter(name=branch_name).delete(permanent=True)
        if child.exists():
            shutil.rmtree(child)
        ln_setup.settings.worktree = False
        ln_setup.settings.dev_dir = previous_dev_dir
        if previous_worktree:
            ln_setup.settings._worktree_path.write_text("true")
