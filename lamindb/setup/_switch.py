from __future__ import annotations

from typing import TYPE_CHECKING

from lamindb_setup import settings

if TYPE_CHECKING:
    from lamindb.models import Branch, Space


def switch(*, branch: str | Branch | None = None, space: str | Space | None = None):
    """Switch to a branch or space, create if not exists."""
    if branch is not None:
        settings.branch = branch
    if space is not None:
        settings.space = space
