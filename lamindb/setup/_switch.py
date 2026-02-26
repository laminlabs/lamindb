from __future__ import annotations

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
    logger.important(f"switched to {target}")
