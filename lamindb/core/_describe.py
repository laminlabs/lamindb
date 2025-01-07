from __future__ import annotations

import datetime
from typing import TYPE_CHECKING

from lamin_utils import logger
from rich.text import Text
from rich.tree import Tree

if TYPE_CHECKING:
    from lamindb.models import Artifact, Collection, Run


def highlight_time(iso: str):
    tz = datetime.datetime.now().astimezone().tzinfo
    res = (
        datetime.datetime.fromisoformat(iso)
        .replace(tzinfo=datetime.timezone.utc)
        .astimezone(tz)
        .strftime("%Y-%m-%d %H:%M:%S")
    )
    return Text(res, style="dim")


# Define consistent column widths
NAME_WIDTH = 25
TYPE_WIDTH = 25
VALUES_WIDTH = 40


def print_rich_tree(tree: Tree, fallback=str):
    from rich.console import Console

    console = Console(force_terminal=True)

    if tree.children:
        try:
            from IPython import get_ipython
            from IPython.core.interactiveshell import InteractiveShell
            from IPython.display import display

            shell = get_ipython()
            if isinstance(shell, InteractiveShell):  # Covers all interactive shells
                display(tree)
                return ""
            else:
                with console.capture() as capture:
                    console.print(tree)
                return capture.get()
        except (ImportError, NameError):
            with console.capture() as capture:
                console.print(tree)
            return capture.get()
    else:
        return fallback


def describe_header(self: Artifact | Collection | Run) -> Tree:
    if hasattr(self, "is_latest") and not self.is_latest:
        logger.warning(
            f"This is not the latest version of the {self.__class__.__name__}."
        )
    if hasattr(self, "_branch_code"):
        if self._branch_code == 0:
            logger.warning("This artifact is hidden.")
        elif self._branch_code == -1:
            logger.warning("This artifact is the trash.")
    # initialize tree
    suffix = self.suffix if hasattr(self, "suffix") and self.suffix else ""
    accessor = self._accessor if hasattr(self, "_accessor") and self._accessor else ""
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
        general.add(
            f".key = '{self.key}'" if self._key_is_virtual else f".key = {self.key}"
        )
    if hasattr(self, "size") and self.size:
        general.add(f".size = {self.size}")
    if hasattr(self, "hash") and self.hash:
        general.add(f".hash = '{self.hash}'")
    if hasattr(self, "n_objects") and self.n_objects:
        general.add(f".n_objects = {self.n_objects}")
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
