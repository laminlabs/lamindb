from __future__ import annotations

import datetime
import re
from typing import TYPE_CHECKING

from django.db import connections
from lamin_utils import colors, logger
from rich.text import Text
from rich.tree import Tree

from .sqlrecord import record_repr

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

    from ..core._context import is_run_from_ipython

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


def describe_artifact_general(
    self: Artifact,
    tree: Tree | None = None,
    foreign_key_data: dict[str, dict[str, int | str]] | None = None,
) -> Tree:
    if tree is None:
        tree = describe_header(self)

    # add general information (order is the same as in API docs)
    general = tree.add(Text("General", style="bold bright_cyan"))

    if self.key:
        general.add(Text.assemble(("key: ", "dim"), (f"{self.key}", "cyan3")))
    if self.description:
        general.add(
            Text.assemble(
                ("description: ", "dim"),
                f"{self.description}",
            )
        )

    # Two column items (short content)
    two_column_items = []

    two_column_items.append(Text.assemble(("uid: ", "dim"), f"{self.uid}"))
    two_column_items.append(Text.assemble(("hash: ", "dim"), f"{self.hash}"))
    two_column_items.append(
        Text.assemble(("size: ", "dim"), f"{format_bytes(self.size)}")
    )
    transform_key = (
        foreign_key_data["run"]["transform_key"]
        if foreign_key_data
        else self.transform.key
        if self.transform
        else None
    )
    two_column_items.append(
        Text.assemble(
            ("transform: ", "dim"),
            (f"{transform_key}", "cyan3"),
        )
    )
    space_name = (
        foreign_key_data["space"]["name"] if foreign_key_data else self.space.name
    )
    two_column_items.append(Text.assemble(("space: ", "dim"), space_name))
    branch_name = (
        foreign_key_data["branch"]["name"] if foreign_key_data else self.space.name
    )
    two_column_items.append(Text.assemble(("branch: ", "dim"), branch_name))
    # actually not name field here, but handle
    created_by_handle = (
        foreign_key_data["created_by"]["name"]
        if foreign_key_data
        else self.created_by.handle
    )
    two_column_items.append(
        Text.assemble(
            ("created_by: ", "dim"),
            (created_by_handle),
        )
    )
    two_column_items.append(
        Text.assemble(("created_at: ", "dim"), highlight_time(str(self.created_at)))
    )
    if self.n_files:
        two_column_items.append(Text.assemble(("n_files: ", "dim"), f"{self.n_files}"))
    if self.n_observations:
        two_column_items.append(
            Text.assemble(("n_observations: ", "dim"), f"{self.n_observations}")
        )
    if self.version:
        two_column_items.append(Text.assemble(("version: ", "dim"), f"{self.version}"))

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

    storage_root = (
        foreign_key_data["storage"]["name"] if foreign_key_data else self.storage.root
    )
    storage_key = self.key if self.key else self.uid
    if not self.key and self.overwrite_versions > 0:
        storage_key = storage_key[:-4]
    general.add(
        Text.assemble(
            ("storage path: ", "dim"),
            (storage_root, "cyan3"),
            f"/{storage_key}",
        )
    )
    return tree


def describe_collection_general(
    self: Collection,
    tree: Tree | None = None,
    foreign_key_data: dict[str, dict[str, int | str]] | None = None,
) -> Tree:
    if tree is None:
        tree = describe_header(self)

    # add general information (order is the same as in API docs)
    general = tree.add(Text("General", style="bold bright_cyan"))

    if self.key:
        general.add(Text.assemble(("key: ", "dim"), (f"{self.key}", "cyan3")))
    if self.description:
        general.add(
            Text.assemble(
                ("description: ", "dim"),
                f"{self.description}",
            )
        )

    # Two column items (short content)
    two_column_items = []

    two_column_items.append(Text.assemble(("uid: ", "dim"), f"{self.uid}"))

    transform_key = (
        foreign_key_data["transform"]["name"]
        if foreign_key_data and "transform" in foreign_key_data
        else self.transform.key
        if self.transform
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
        else self.space.name
        if self.space
        else None
    )
    if space_name:
        two_column_items.append(Text.assemble(("space: ", "dim"), space_name))

    branch_name = (
        foreign_key_data["branch"]["name"]
        if foreign_key_data and "branch" in foreign_key_data
        else self.branch.name
        if self.branch
        else None
    )
    if branch_name:
        two_column_items.append(Text.assemble(("branch: ", "dim"), branch_name))

    created_by_handle = (
        foreign_key_data["created_by"]["name"]
        if foreign_key_data and "created_by" in foreign_key_data
        else self.created_by.handle
        if self.created_by
        else None
    )
    if created_by_handle:
        two_column_items.append(
            Text.assemble(
                ("created_by: ", "dim"),
                (created_by_handle),
            )
        )

    if self.created_at:
        two_column_items.append(
            Text.assemble(("created_at: ", "dim"), highlight_time(str(self.created_at)))
        )

    if self.version:
        two_column_items.append(Text.assemble(("version: ", "dim"), f"{self.version}"))

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
            padding_needed = max(0, 45 - len(left_plain_text))
            padding = " " * padding_needed

            general.add(Text.assemble(left_item, padding, right_item))
        else:
            # Single item (odd number)
            general.add(two_column_items[i])

    return tree


def describe_run_general(
    self: Run,
    tree: Tree | None = None,
    foreign_key_data: dict[str, dict[str, int | str]] | None = None,
) -> Tree:
    if tree is None:
        tree = describe_header(self)

    # add general information (order is the same as in API docs)
    general = tree.add(Text("General", style="bold bright_cyan"))

    if self.name:
        general.add(Text.assemble(("name: ", "dim"), (f"{self.name}", "cyan3")))

    # Two column items (short content)
    two_column_items = []

    two_column_items.append(Text.assemble(("uid: ", "dim"), f"{self.uid}"))

    transform_key = (
        foreign_key_data["transform"]["name"]
        if foreign_key_data and "transform" in foreign_key_data
        else self.transform.key
        if self.transform
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
        else self.space.name
        if self.space
        else None
    )
    if space_name:
        two_column_items.append(Text.assemble(("space: ", "dim"), space_name))

    branch_name = (
        foreign_key_data["branch"]["name"]
        if foreign_key_data and "branch" in foreign_key_data
        else self.branch.name
        if self.branch
        else None
    )
    if branch_name:
        two_column_items.append(Text.assemble(("branch: ", "dim"), branch_name))

    created_by_handle = (
        foreign_key_data["created_by"]["name"]
        if foreign_key_data and "created_by" in foreign_key_data
        else self.created_by.handle
        if self.created_by
        else None
    )
    if created_by_handle:
        two_column_items.append(
            Text.assemble(
                ("created_by: ", "dim"),
                (created_by_handle),
            )
        )

    if self.created_at:
        two_column_items.append(
            Text.assemble(("created_at: ", "dim"), highlight_time(str(self.created_at)))
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
            padding_needed = max(0, 45 - len(left_plain_text))
            padding = " " * padding_needed

            general.add(Text.assemble(left_item, padding, right_item))
        else:
            # Single item (odd number)
            general.add(two_column_items[i])
    if self.params:
        params = tree.add(Text("Params", style="bold bright_cyan"))
        for key, value in self.params.items():
            params.add(f"{key}: {value}")
    return tree


def describe_postgres(self):
    from ._django import get_artifact_or_run_with_related, get_collection_with_related
    from ._feature_manager import describe_features

    model_name = self.__class__.__name__
    msg = f"{colors.green(model_name)}{record_repr(self, include_foreign_keys=False).lstrip(model_name)}\n"
    if self._state.db is not None and self._state.db != "default":
        msg += f"  {colors.italic('Database instance')}\n"
        msg += f"    slug: {self._state.db}\n"

    if model_name in {"Artifact", "Run"}:
        result = get_artifact_or_run_with_related(
            self,
            include_feature_link=True,
            include_fk=True,
            include_m2m=True,
            include_schema=True,
        )
        related_data = result.get("related_data", {})
        if model_name == "Artifact":
            tree = describe_artifact_general(self, foreign_key_data=related_data["fk"])
        else:
            tree = describe_run_general(self, foreign_key_data=related_data["fk"])
        return describe_features(
            self,
            tree=tree,
            related_data=related_data,
            with_labels=True,
        )
    elif model_name == "Collection":
        result = get_collection_with_related(self, include_fk=True)
        related_data = result.get("related_data", {})
        tree = describe_collection_general(
            self, foreign_key_data=related_data.get("fk", {})
        )
        return tree
    else:
        tree = describe_header(self)
        return tree


def describe_sqlite(self, print_types: bool = False):
    from ._feature_manager import describe_features

    model_name = self.__class__.__name__
    msg = f"{colors.green(model_name)}{record_repr(self, include_foreign_keys=False).lstrip(model_name)}\n"
    if self._state.db is not None and self._state.db != "default":
        msg += f"  {colors.italic('Database instance')}\n"
        msg += f"    slug: {self._state.db}\n"

    fields = self._meta.fields
    direct_fields = []
    foreign_key_fields = []
    for f in fields:
        if f.is_relation:
            foreign_key_fields.append(f.name)
        else:
            direct_fields.append(f.name)
    if not self._state.adding:
        # prefetch foreign key relationships
        self = (
            self.__class__.objects.using(self._state.db)
            .select_related(*foreign_key_fields)
            .get(id=self.id)
        )
        # prefetch m-2-m relationships
        many_to_many_fields = []
        if model_name in {"Artifact", "Collection"}:
            many_to_many_fields.append("input_of_runs")
        if model_name == "Artifact":
            many_to_many_fields.append("feature_sets")
        self = (
            self.__class__.objects.using(self._state.db)
            .prefetch_related(*many_to_many_fields)
            .get(id=self.id)
        )
    if model_name in {"Artifact", "Run"}:
        if model_name == "Artifact":
            tree = describe_artifact_general(self)
        else:
            tree = describe_run_general(self)
        return describe_features(
            self,
            tree=tree,
            with_labels=True,
        )
    elif model_name == "Collection":
        tree = describe_collection_general(self)
        return tree
    else:
        tree = describe_header(self)
        return tree


def describe_postgres_sqlite(self, return_str: bool = False) -> str | None:
    from ._describe import format_rich_tree

    if not self._state.adding and connections[self._state.db].vendor == "postgresql":
        tree = describe_postgres(self)
    else:
        tree = describe_sqlite(self)

    return format_rich_tree(tree, return_str=return_str)
