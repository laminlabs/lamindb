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
    return Text(res, style="dim")


# Define consistent column widths
NAME_WIDTH = 25
TYPE_WIDTH = 25
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
    suffix_accessor = (
        f"{suffix}/{accessor}" if suffix and accessor else suffix or accessor or ""
    )

    tree = Tree(
        Text.assemble(
            (self.__class__.__name__, "bold"), (f" {suffix_accessor}", "bold dim")
        ),
        guide_style="dim",  # dim the connecting lines
    )
    return tree


def describe_general(self: Artifact | Collection, tree: Tree | None = None) -> Tree:
    if tree is None:
        tree = describe_header(self)

    # add general information (order is the same as in API docs)
    general = tree.add(Text("General", style="bold bright_cyan"))
    general.add(f".uid = '{self.uid}'")
    if hasattr(self, "key") and self.key:
        general.add(f".key = '{self.key}'")
    if hasattr(self, "size") and self.size:
        general.add(f".size = {self.size}")
    if hasattr(self, "hash") and self.hash:
        general.add(f".hash = '{self.hash}'")
    if hasattr(self, "n_files") and self.n_files:
        general.add(f".n_files = {self.n_files}")
    if hasattr(self, "n_observations") and self.n_observations:
        general.add(Text(f".n_observations = {self.n_observations}"))
    if hasattr(self, "version") and self.version:
        general.add(Text(f".version = '{self.version}'"))

    if hasattr(self, "storage"):
        storage_root = self.storage.root
        # general.add(f".storage = {storage_root}")
        general.add(
            Text.assemble(
                ".path = ",
                (storage_root, "dim"),
                f"{str(self.path).removeprefix(storage_root)}",
            )
        )
    if hasattr(self, "created_by") and self.created_by:
        general.add(
            Text.assemble(
                ".created_by = ",
                (
                    self.created_by.handle
                    if self.created_by.name is None
                    else f"{self.created_by.handle} ({self.created_by.name})"
                ),
            )
        )
    if hasattr(self, "created_at") and self.created_at:
        general.add(
            Text.assemble(".created_at = ", highlight_time(str(self.created_at)))
        )
    if hasattr(self, "transform") and self.transform:
        general.add(
            Text(
                f".transform = '{self.transform.description}'",
                style="cyan3",
            )
        )

    return tree
