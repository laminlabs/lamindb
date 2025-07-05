from __future__ import annotations

import datetime
import re
from typing import TYPE_CHECKING

from lamin_utils import logger
from rich.text import Text
from rich.tree import Tree

from ..core._context import is_run_from_ipython

if TYPE_CHECKING:
    from lamindb.models import Artifact, Collection, Run


def highlight_time(iso: str):
    try:
        tz = datetime.datetime.now().astimezone().tzinfo
        res = (
            datetime.datetime.fromisoformat(iso)
            .replace(tzinfo=datetime.timezone.utc)
            .astimezone(tz)
            .strftime("%Y-%m-%d %H:%M:%S")
        )
    except ValueError:
        # raises ValueError: Invalid isoformat string: '<django.db.models.expressions.DatabaseDefault object at 0x1128ac440>'
        # but can't be caught with `isinstance(iso, DatabaseDefault)` for unkown reasons
        return Text("timestamp of unsaved record not available", style="dim")
    return Text(res)


# Define consistent column widths
NAME_WIDTH = 30
TYPE_WIDTH = 35  # types can get long, e.g. cat[Record[Treatment]]
VALUES_WIDTH = 40


def format_rich_tree(
    tree: Tree, fallback: str = "", return_str: bool = False, strip_ansi: bool = True
) -> str | None:
    from rich.console import Console

    # If tree has no children, return fallback
    if not tree.children:
        return fallback

    console = Console(force_terminal=True)
    printed = False

    if return_str:
        from io import StringIO

        string_io = StringIO()
        str_console = Console(file=string_io, force_terminal=True)
        str_console.print(tree)
        result = string_io.getvalue()
        if strip_ansi:
            ansi_escape = re.compile(r"\x1b(?:\[[0-9;]*[a-zA-Z]|\(B)")
            result = ansi_escape.sub("", result)
        # rstrip trailing whitespace on every line
        result = "\n".join(line.rstrip() for line in result.splitlines())
        return result

    try:
        if not is_run_from_ipython:
            from IPython import get_ipython
            from IPython.core.interactiveshell import InteractiveShell
            from IPython.display import display

            shell = get_ipython()
            if isinstance(shell, InteractiveShell):
                display(tree)
                printed = True
                return None
    except (NameError, ImportError):
        pass

    if not printed:
        # be careful to test this on a terminal
        console = Console(force_terminal=True)
        console.print(tree)

    return None


def describe_header(self: Artifact | Collection | Run) -> Tree:
    if hasattr(self, "is_latest") and not self.is_latest:
        logger.warning(
            f"This is not the latest version of the {self.__class__.__name__}."
        )
    if hasattr(self, "branch_id"):
        if self.branch_id == 0:  # type: ignore
            logger.warning("This artifact is archived.")
        elif self.branch_id == -1:  # type: ignore
            logger.warning("This artifact is in the trash.")
    # initialize tree
    suffix = self.suffix if hasattr(self, "suffix") and self.suffix else ""
    accessor = self.otype if hasattr(self, "otype") and self.otype else ""
    kind = f" · {self.kind}" if hasattr(self, "kind") and self.kind else ""
    suffix_accessor = (
        f"{suffix} · {accessor}{kind}"
        if suffix and accessor
        else suffix or accessor or ""
    )

    tree = Tree(
        Text.assemble(
            (self.__class__.__name__, "bold"), (f" {suffix_accessor}", "bold dim")
        ),
        guide_style="dim",  # dim the connecting lines
    )
    return tree


def format_bytes(bytes_value):
    """Convert bytes to human readable format."""
    if bytes_value < 1024:
        return f"{bytes_value} B"
    elif bytes_value < 1024**2:
        return f"{bytes_value / 1024:.1f} KB"
    elif bytes_value < 1024**3:
        return f"{bytes_value / (1024**2):.1f} MB"
    elif bytes_value < 1024**4:
        return f"{bytes_value / (1024**3):.1f} GB"
    else:
        return f"{bytes_value / (1024**4):.1f} TB"


def describe_general(self: Artifact | Collection, tree: Tree | None = None) -> Tree:
    if tree is None:
        tree = describe_header(self)

    # add general information (order is the same as in API docs)
    general = tree.add(Text("General", style="bold bright_cyan"))

    # Two column items (short content)
    two_column_items = []

    two_column_items.append(Text.assemble(("uid: ", "dim"), f"{self.uid}"))
    if hasattr(self, "hash") and self.hash:
        two_column_items.append(Text.assemble(("hash: ", "dim"), f"{self.hash}"))
    if hasattr(self, "size") and self.size:
        two_column_items.append(
            Text.assemble(("size: ", "dim"), f"{format_bytes(self.size)}")
        )
    if hasattr(self, "n_files") and self.n_files:
        two_column_items.append(Text.assemble(("n_files: ", "dim"), f"{self.n_files}"))
    if hasattr(self, "n_observations") and self.n_observations:
        two_column_items.append(
            Text.assemble(("n_observations: ", "dim"), f"{self.n_observations}")
        )
    if hasattr(self, "version") and self.version:
        two_column_items.append(Text.assemble(("version: ", "dim"), f"{self.version}"))
    if hasattr(self, "space"):
        two_column_items.append(Text.assemble(("space: ", "dim"), f"{self.space.name}"))
    if hasattr(self, "branch"):
        two_column_items.append(
            Text.assemble(("branch: ", "dim"), f"{self.branch.name}")
        )
    if hasattr(self, "created_at") and self.created_at:
        two_column_items.append(
            Text.assemble(("created_at: ", "dim"), highlight_time(str(self.created_at)))
        )
    if hasattr(self, "created_by") and self.created_by:
        two_column_items.append(
            Text.assemble(
                ("created_by: ", "dim"),
                (
                    self.created_by.handle
                    if self.created_by.name is None
                    else f"{self.created_by.handle} ({self.created_by.name})"
                ),
            )
        )

    # Add two-column items in pairs
    for i in range(0, len(two_column_items), 2):
        if i + 1 < len(two_column_items):
            # Two items side by side
            left_item = two_column_items[i]
            right_item = two_column_items[i + 1]

            # Create padded version by calculating the plain text length
            left_plain_text = (
                left_item.plain if hasattr(left_item, "plain") else str(left_item)
            )
            padding_needed = max(0, 35 - len(left_plain_text))
            padding = " " * padding_needed

            general.add(Text.assemble(left_item, padding, right_item))
        else:
            # Single item (odd number)
            general.add(two_column_items[i])

    # Single column items (long content)
    if hasattr(self, "key") and self.key:
        general.add(Text.assemble(("key: ", "dim"), f"{self.key}"))
    if hasattr(self, "storage"):
        storage_root = self.storage.root
        general.add(
            Text.assemble(
                ("storage location / path: ", "dim"),
                (storage_root, "cyan3"),
                f"{str(self.path).removeprefix(storage_root)}",
            )
        )
    if hasattr(self, "description") and self.description is not None:
        general.add(
            Text.assemble(
                ("description: ", "dim"),
                f"{self.description}",
            )
        )
    if hasattr(self, "transform") and self.transform is not None:
        general.add(
            Text.assemble(
                ("transform: ", "dim"),
                f"{self.transform.key}",
            )
        )

    return tree
