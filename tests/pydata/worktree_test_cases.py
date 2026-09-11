from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import lamindb as ln
import lamindb_setup as ln_setup
from lamindb_setup.core._settings_store import local_current_branch_file
from lamindb_setup.errors import WorktreePathError


def transform_dev_dir_key(root: Path) -> None:
    path = root / "pipelines" / f"wf-{time.time_ns()}.py"
    path.parent.mkdir(parents=True)
    path.write_text("print('hello')\n")
    ln_setup.settings.dev_dir = root
    assert ln.Transform.from_path(path).key == f"pipelines/{path.name}"


def transform_relative_dev_dir_key(root: Path) -> None:
    path = root / "pipelines" / f"wf-{time.time_ns()}.py"
    path.parent.mkdir(parents=True)
    path.write_text("print('hello')\n")
    ln_setup.settings.dev_dir = root
    os.chdir(root)
    assert ln.Transform.from_path(Path("pipelines") / path.name).key == (
        f"pipelines/{path.name}"
    )


def transform_worktree_key(root: Path) -> None:
    root.mkdir(parents=True)
    ln_setup.settings.dev_dir = root
    ln_setup.settings.worktree = True
    child = root / "feature-a"
    path = child / "pipelines" / f"wf-{time.time_ns()}.py"
    path.parent.mkdir(parents=True)
    path.write_text("print('hello from worktree')\n")
    os.chdir(child)
    assert ln.Transform.from_path(path).key == f"pipelines/{path.name}"


def transform_outside_worktree_errors(root: Path) -> None:
    root.mkdir(parents=True)
    ln_setup.settings.dev_dir = root
    ln_setup.settings.worktree = True
    child = root / "feature-a"
    path = child / "pipelines" / f"wf-{time.time_ns()}.py"
    path.parent.mkdir(parents=True)
    path.write_text("print('outside child should fail')\n")
    os.chdir(root)
    try:
        ln.Transform.from_path(path)
    except WorktreePathError as error:
        assert "inside a child directory" in str(error)
    else:
        raise AssertionError("a path outside a worktree child should be rejected")


def switch_creates_worktree(root: Path) -> None:
    root.mkdir(parents=True)
    ln_setup.settings.dev_dir = root
    ln_setup.settings.worktree = True
    branch_name = f"wt-{time.time_ns()}"
    os.chdir(root)
    try:
        ln.setup.switch(branch_name, create=True)
        marker = local_current_branch_file(root / branch_name)
        assert marker.read_text().split("\n")[1] == branch_name
    finally:
        ln.Branch.filter(name=branch_name).delete(permanent=True)


if __name__ == "__main__":
    case_name, path = sys.argv[1:]
    globals()[case_name](Path(path))
