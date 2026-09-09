from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

from lamin_utils import logger
from lamindb_setup import settings

if TYPE_CHECKING:
    from lamindb.models import Branch


def switch(target: str | Branch, *, space: bool = False, create: bool = False):
    """Switch to a branch or space, create if not exists.

    Args:
        target: Branch target or space target to switch to.
        space: If True, switch space; otherwise switch branch.
        create: If True and switching branch, create the branch if it does not exist.
    """
    if space:
        settings.space = target
    else:
        if settings.worktree:
            if isinstance(target, str) and "/" in target:
                raise ValueError(
                    "Branch names containing '/' are not supported in worktree mode."
                )
            dev_dir = settings.dev_dir
            if dev_dir is not None:
                dev_dir = dev_dir.resolve()
                cwd = Path.cwd().resolve()
                if create and cwd == dev_dir and isinstance(target, str):
                    child_dir = dev_dir / target
                    if child_dir.exists() and not child_dir.is_dir():
                        raise ValueError(
                            f"Cannot create worktree directory '{child_dir}': path exists and is not a directory."
                        )
                    child_dir.mkdir(parents=True, exist_ok=True)
                    original_cwd = cwd
                    try:
                        os.chdir(child_dir)
                        switch(target, space=False, create=True)
                    finally:
                        os.chdir(original_cwd)
                    return
            settings._resolve_active_worktree_root(raise_on_error=True)

        is_worktree_bootstrap = (
            create
            and settings.worktree
            and isinstance(target, str)
            and settings.dev_dir is not None
            and Path.cwd().resolve().parent == settings.dev_dir.resolve()
            and Path.cwd().resolve().name == target
            and not settings._branch_path.exists()
        )
        if create:
            from lamindb import Branch, Q
            from lamindb.errors import BranchAlreadyExists

            # Consistent with git switch -c: error if branch already exists.
            existing = Branch.filter(Q(name=target) | Q(uid=target)).one_or_none()
            if existing is not None:
                raise BranchAlreadyExists(
                    f"Branch '{target}' already exists. Omit -c/--create to switch to it."
                )
            Branch(name=target).save()
            logger.important(f"created branch: {target}")
        settings.branch = target
    if is_worktree_bootstrap:
        logger.important_hint(f"to switch, cd into {target}")
    else:
        logger.important(f"switched to {target}")
