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


def describe_dataset(
    obj: Artifact | Collection,
    tree: Tree | None = None,
    foreign_key_data: dict[str, dict[str, int | str]] | None = None,
) -> Tree:
    from lamindb.models import Artifact

    if tree is None:
        tree = describe_header(obj)

    general = tree.add(Text("General", style="bold bright_cyan"))

    if hasattr(obj, "key") and obj.key:
        general.add(Text.assemble(("key: ", "dim"), (f"{obj.key}", "cyan3")))
    if hasattr(obj, "description") and obj.description:
        general.add(
            Text.assemble(
                ("description: ", "dim"),
                f"{obj.description}",
            )
        )

    two_column_items = []
    two_column_items.append(Text.assemble(("uid: ", "dim"), f"{obj.uid}"))

    if isinstance(obj, Artifact):
        two_column_items.append(Text.assemble(("hash: ", "dim"), f"{obj.hash}"))
        two_column_items.append(
            Text.assemble(("size: ", "dim"), f"{format_bytes(obj.size)}")
        )

    transform_key = (
        foreign_key_data["run"]["transform_key"]
        if foreign_key_data and isinstance(obj, Artifact)
        else foreign_key_data["transform"]["name"]
        if foreign_key_data and "transform" in foreign_key_data
        else obj.transform.key
        if hasattr(obj, "transform") and obj.transform and isinstance(obj, Artifact)
        else obj.transform.name
        if hasattr(obj, "transform") and obj.transform
        else None
    )
    if transform_key:
        two_column_items.append(
            Text.assemble(
                ("transform: ", "dim"),
                (f"{transform_key}", "cyan3"),
            )
        )

    space_name = (
        foreign_key_data["space"]["name"]
        if foreign_key_data and "space" in foreign_key_data
        else obj.space.name
        if hasattr(obj, "space") and obj.space
        else None
    )
    if space_name:
        two_column_items.append(Text.assemble(("space: ", "dim"), space_name))

    branch_name = (
        foreign_key_data["branch"]["name"]
        if foreign_key_data and "branch" in foreign_key_data
        else getattr(
            obj.space if isinstance(obj, Artifact) else obj.branch, "name", None
        )
        if hasattr(obj, "branch" if hasattr(obj, "branch") else "space")
        else None
    )
    if branch_name:
        two_column_items.append(Text.assemble(("branch: ", "dim"), branch_name))

    created_by_handle = (
        foreign_key_data["created_by"]["name"]
        if foreign_key_data and "created_by" in foreign_key_data
        else foreign_key_data["branch"]["name"]
        if foreign_key_data and isinstance(obj, Artifact)
        else obj.created_by.handle
        if hasattr(obj, "created_by") and obj.created_by
        else None
    )
    if created_by_handle:
        two_column_items.append(
            Text.assemble(
                ("created_by: ", "dim"),
                (created_by_handle),
            )
        )

    if hasattr(obj, "created_at") and obj.created_at:
        two_column_items.append(
            Text.assemble(("created_at: ", "dim"), highlight_time(str(obj.created_at)))
        )

    if isinstance(obj, Artifact):
        if obj.n_files:
            two_column_items.append(
                Text.assemble(("n_files: ", "dim"), f"{obj.n_files}")
            )
        if obj.n_observations:
            two_column_items.append(
                Text.assemble(("n_observations: ", "dim"), f"{obj.n_observations}")
            )

    if hasattr(obj, "version") and obj.version:
        two_column_items.append(Text.assemble(("version: ", "dim"), f"{obj.version}"))

    for i in range(0, len(two_column_items), 2):
        if i + 1 < len(two_column_items):
            # Two items side by side
            left_item = two_column_items[i]
            right_item = two_column_items[i + 1]

            # Create padded version by calculating the plain text length
            left_plain_text = (
                left_item.plain if hasattr(left_item, "plain") else str(left_item)
            )
            padding_needed = max(0, 45 - len(left_plain_text))
            padding = " " * padding_needed

            general.add(Text.assemble(left_item, padding, right_item))
        else:
            # Single item (odd number)
            general.add(two_column_items[i])

    # Should appear at the end so we need another if statement
    if isinstance(obj, Artifact):
        storage_root = (
            foreign_key_data["storage"]["name"]
            if foreign_key_data
            else obj.storage.root
        )
        storage_key = obj.key if obj.key else obj.uid
        if not obj.key and obj.overwrite_versions > 0:
            storage_key = storage_key[:-4]
        general.add(
            Text.assemble(
                ("storage path: ", "dim"),
                (storage_root, "cyan3"),
                f"/{storage_key}",
            )
        )

    return tree
