from __future__ import annotations

import re
from types import SimpleNamespace
from typing import TYPE_CHECKING

from django.db import connections
from lamin_utils import colors, logger
from rich.table import Column, Table
from rich.text import Text
from rich.tree import Tree

from lamindb.models import Run

from .sqlrecord import SQLRecord, format_field_value, record_repr

if TYPE_CHECKING:
    from lamindb.models import Artifact, Collection, Schema, Transform

    from .run import TracksRun


# Define consistent column widths for use in other modules
NAME_WIDTH = 30
TYPE_WIDTH = 35  # types can get long, e.g. cat[Record[Treatment]]
VALUES_WIDTH = 40


def strip_ansi_from_string(text: str) -> str:
    """Remove ANSI escape sequences from a string."""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


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
            result = strip_ansi_from_string(result)
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


def format_run_title(
    record: Run | SimpleNamespace | None,
    transform_key: str | None = None,
    dim: bool = False,
) -> Text:
    if record is None:
        return Text("")
    display_name = (
        Text(record.name, style="cyan3")
        if record.name is not None
        else Text(record.uid[:7], style="cyan3")
    )
    if transform_key is None:
        transform_key = record.transform.key
    title = Text.assemble(
        display_name,
        (" (", "dim"),
        (transform_key, "cyan3"),
        (")", "dim"),
    )
    return title


def format_title_with_version(
    record: Artifact | Collection | Transform | SimpleNamespace,
) -> Text:
    title_str = record.key if record.key is not None else ""
    title = Text.assemble(
        (title_str, "cyan3"),
        (f" ({record.version if record.version else record.uid[-4:]})", "dim"),
        Text.assemble(("\n|   description: ", "dim"), record.description)
        if record.description
        else Text(""),
    )
    return title


def describe_header(record: Artifact | Collection | Run | Transform) -> Tree:
    if hasattr(record, "is_latest") and not record.is_latest:
        logger.warning(
            f"This is not the latest version of the {record.__class__.__name__}."
        )
    if record.branch_id == 0:  # type: ignore
        logger.warning("This artifact is archived.")
    elif record.branch_id == -1:  # type: ignore
        logger.warning("This artifact is in the trash.")
    if isinstance(record, Run):
        title = format_run_title(record, dim=True)  # dim makes the uid grey
    else:
        title = format_title_with_version(record)
    tree = Tree(
        Text.assemble(
            (f"{record.__class__.__name__}: ", "bold"),
            title,
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


def append_uid_run(record: TracksRun, two_column_items: list, fk_data=None) -> None:
    if fk_data and "run" in fk_data and fk_data["run"] and fk_data["run"]["id"]:
        run, transform_key = (
            SimpleNamespace(**fk_data["run"]),
            fk_data["run"]["transform_key"],
        )
    elif record.run is not None:
        run, transform_key = record.run, record.run.transform.key
    else:
        run, transform_key = None, None
    text_uid = Text.assemble(("uid: ", "dim"), f"{record.uid}")
    text_run = Text.assemble(
        ("run: ", "dim"), format_run_title(run, transform_key=transform_key)
    )
    two_column_items.append(text_uid)
    two_column_items.append(text_run)


def append_branch_space_created_at_created_by(
    record: SQLRecord, two_column_items, fk_data=None
):
    # branch
    branch_name = fk_data["branch"]["name"] if fk_data else record.branch.name
    two_column_items.append(Text.assemble(("branch: ", "dim"), branch_name))
    # space
    space_name = fk_data["space"]["name"] if fk_data else record.space.name
    two_column_items.append(Text.assemble(("space: ", "dim"), space_name))
    # created_at
    two_column_items.append(
        Text.assemble(("created_at: ", "dim"), format_field_value(record.created_at))
    )
    # created_by / "name" in fk_data holds handle, is display name
    created_by_handle = (
        fk_data["created_by"]["name"] if fk_data else record.created_by.handle
    )
    two_column_items.append(Text.assemble(("created_by: ", "dim"), created_by_handle))


def add_two_column_items_to_tree(tree: Tree, two_column_items: list) -> None:
    table = Table(
        Column("", no_wrap=True),
        Column("", no_wrap=True),
        show_header=False,
        box=None,
        pad_edge=False,
    )
    for i in range(0, len(two_column_items), 2):
        if i + 1 < len(two_column_items):
            left_item = two_column_items[i]
            right_item = two_column_items[i + 1]
            table.add_row(left_item, right_item)
        else:
            table.add_row(two_column_items[i], "")
    tree.add(table)


def describe_artifact(
    record: Artifact,
    related_data: dict | None = None,
) -> Tree:
    from ._feature_manager import describe_features
    from ._label_manager import describe_labels

    if related_data is not None:
        fk_data = related_data.get("fk", {})
    else:
        fk_data = {}
    tree = describe_header(record)
    dataset_features_tree, external_features_tree = describe_features(
        record,
        related_data=related_data,
    )
    labels_tree = describe_labels(record, related_data=related_data)
    two_column_items = []  # type: ignore
    append_uid_run(record, two_column_items, fk_data)
    if record.kind or record.otype:
        two_column_items.append(Text.assemble(("kind: ", "dim"), f"{record.kind}"))
        two_column_items.append(Text.assemble(("otype: ", "dim"), f"{record.otype}"))
    two_column_items.append(Text.assemble(("hash: ", "dim"), f"{record.hash}"))
    two_column_items.append(
        Text.assemble(("size: ", "dim"), f"{format_bytes(record.size)}")
    )
    append_branch_space_created_at_created_by(record, two_column_items, fk_data)
    if record.n_files:
        two_column_items.append(
            Text.assemble(("n_files: ", "dim"), f"{record.n_files}")
        )
    if record.n_observations:
        two_column_items.append(
            Text.assemble(("n_observations: ", "dim"), f"{record.n_observations}")
        )
    add_two_column_items_to_tree(tree, two_column_items)
    storage_root = fk_data["storage"]["name"] if fk_data else record.storage.root
    storage_key = (
        record.key
        if not record._key_is_virtual
        else record._real_key
        if record._real_key
        else f".lamindb/{record.uid}"
    )
    if record.uid in storage_key:
        if record.overwrite_versions:
            storage_key = storage_key[:-4]
        storage_key = f"{storage_key}{record.suffix}"
    tree.add(
        Text.assemble(
            ("storage/path: ", "dim"),
            (storage_root, "cyan3"),
            ("/", "dim"),
            storage_key,
        )
    )
    if dataset_features_tree:
        tree.add(dataset_features_tree)
    if external_features_tree:
        tree.add(external_features_tree)
    if labels_tree:
        tree.add(labels_tree)
    return tree


def describe_collection(
    record: Collection,
    related_data: dict | None = None,
) -> Tree:
    tree = describe_header(record)
    if related_data is not None:
        fk_data = related_data.get("fk", {})
    else:
        fk_data = {}
    two_column_items = []  # type: ignore
    append_uid_run(record, two_column_items, fk_data)
    append_branch_space_created_at_created_by(record, two_column_items, fk_data)
    add_two_column_items_to_tree(tree, two_column_items)
    return tree


def display_text(
    text: str, title: str, tree: Tree, max_lines: int = 30, uid: str = ""
) -> None:
    # Split the code into lines and add dim vertical bars
    lines = text.split("\n")
    end_parts = [("\n│ …", "grey30")] if len(lines) > max_lines else []
    parts = [(title + ": ", "purple")]
    parts.append((uid, ""))
    max_length = 80
    for line in lines[:max_lines]:
        parts.append(("\n│ ", "dim"))
        parts.append((line[:max_length], "grey30"))
        if len(line) > max_length:
            parts.append((" …", "grey30"))
    parts.extend(end_parts)
    tree.add(Text.assemble(*parts))


def describe_run(
    record: Run,
    related_data: dict | None = None,
) -> Tree:
    from ._feature_manager import describe_features

    tree = describe_header(record)
    if related_data is not None:
        fk_data = related_data.get("fk", {})
    else:
        fk_data = {}
    _, features_tree = describe_features(
        record,
        related_data=related_data,
    )
    two_column_items = []  # type: ignore
    two_column_items.append(Text.assemble(("uid: ", "dim"), f"{record.uid}"))
    if fk_data and "transform" in fk_data:
        transform = SimpleNamespace(**fk_data["transform"], description="")
    else:
        transform = record.transform
    two_column_items.append(
        Text.assemble(
            ("transform: ", "dim"),
            format_title_with_version(transform),
        )
    )
    two_column_items.append(
        Text.assemble(
            ("started_at: ", "dim"), format_field_value(record.started_at, none="")
        )
    )
    two_column_items.append(
        Text.assemble(
            ("finished_at: ", "dim"), format_field_value(record.finished_at, none="")
        )
    )
    two_column_items.append(Text.assemble(("status: ", "dim"), record.status))
    two_column_items.append(
        Text.assemble(("reference: ", "dim"), record.reference)
        if record.reference
        else Text("")
    )
    append_branch_space_created_at_created_by(record, two_column_items, fk_data)
    add_two_column_items_to_tree(tree, two_column_items)
    if record.cli_args:
        display_text(
            record.cli_args.strip(),
            "cli_args",
            tree,
            max_lines=4,
        )
    if record.report_id:
        report = record.report.load(is_run_input=False)
        if report:
            display_text(
                strip_ansi_from_string(report.strip()),
                "report",
                tree,
                max_lines=4,
                uid=record.report.uid[:7],
            )
    if record.environment_id:
        display_text(
            record.environment.load(is_run_input=False).strip(),
            "environment",
            tree,
            max_lines=4,
            uid=record.environment.uid[:7],
        )
    if record.params:
        params = tree.add(Text("Params", style="bold dark_orange"))
        for key, value in record.params.items():
            params.add(f"{key}: {value}")
    if features_tree:
        tree.add(features_tree)
    return tree


def describe_transform(
    record: Transform,
    related_data: dict | None = None,
) -> Tree:
    tree = describe_header(record)
    if related_data is not None:
        fk_data = related_data.get("fk", {})
    else:
        fk_data = {}
    two_column_items = []  # type: ignore
    two_column_items.append(Text.assemble(("uid: ", "dim"), f"{record.uid}"))
    two_column_items.append(
        Text.assemble(("reference: ", "dim"), record.reference)
        if record.reference
        else Text("")
    )
    two_column_items.append(Text.assemble(("hash: ", "dim"), f"{record.hash}"))
    two_column_items.append(Text.assemble(("type: ", "dim"), f"{record.type}"))
    append_branch_space_created_at_created_by(record, two_column_items, fk_data)
    add_two_column_items_to_tree(tree, two_column_items)
    if record.source_code:
        display_text(record.source_code.strip(), "source_code", tree)
    return tree


def describe_schema(record: Schema, slot: str | None = None) -> Tree:
    from ._feature_manager import strip_cat

    if record.type:
        prefix = f" {record.type.name} · "
    else:
        prefix = " "
    if record.name:
        name = record.name
    else:
        name = "unnamed"
    header = "Schema:" if slot is None else f"{slot}:"
    description = (
        Text.assemble(("\n|   description: ", "dim"), record.description)
        if record.description
        else Text("")
    )
    tree = Tree(
        Text.assemble(
            (header, "bold"), (f"{prefix}", "dim"), (f"{name}", "cyan3"), description
        ),
        guide_style="dim",
    )
    two_column_items = []  # type: ignore
    append_uid_run(record, two_column_items)
    two_column_items.append(Text.assemble(("itype: ", "dim"), f"{record.itype}"))
    two_column_items.append(Text.assemble(("otype: ", "dim"), f"{record.otype}"))
    two_column_items.append(Text.assemble(("hash: ", "dim"), f"{record.hash}"))
    two_column_items.append(
        Text.assemble(("ordered_set: ", "dim"), f"{record.ordered_set}")
    )
    two_column_items.append(
        Text.assemble(("maximal_set: ", "dim"), f"{record.maximal_set}")
    )
    two_column_items.append(
        Text.assemble(("minimal_set: ", "dim"), f"{record.minimal_set}")
    )
    append_branch_space_created_at_created_by(record, two_column_items)
    add_two_column_items_to_tree(tree, two_column_items)

    # Add features section
    members_count = record.n
    members_count_display = f" ({members_count})" if members_count > 0 else ""
    if members_count > 0 or (record.dtype and record.itype is not None):
        features = tree.add(
            Text.assemble(
                (
                    "Features" if record.itype == "Feature" else record.itype,
                    "bold bright_magenta",
                ),
                (members_count_display, "bold dim"),
            )
        )
        if members_count > 0:
            feature_table = Table(
                show_header=True, header_style="dim", box=None, pad_edge=False
            )

            feature_table.add_column("name", style="", no_wrap=True)
            feature_table.add_column("dtype", style="", no_wrap=True)
            feature_table.add_column("optional", style="", no_wrap=True)
            feature_table.add_column("nullable", style="", no_wrap=True)
            feature_table.add_column("coerce_dtype", style="", no_wrap=True)
            feature_table.add_column("default_value", style="", no_wrap=True)

            optionals = record.optionals.get()
            for member in record.members:
                feature_table.add_row(
                    member.name,
                    Text(strip_cat(member.dtype)),
                    "✓" if optionals.filter(uid=member.uid).exists() else "✗",
                    "✓" if member.nullable else "✗",
                    "✓" if record.coerce_dtype or member.coerce_dtype else "✗",
                    str(member.default_value) if member.default_value else "unset",
                )

            features.add(feature_table)
        elif record.dtype:
            features.add(Text.assemble(("dtype: ", "dim"), f"{record.dtype}"))

    return tree


def describe_postgres(record):
    from ._django import get_artifact_or_run_with_related, get_collection_with_related

    model_name = record.__class__.__name__
    msg = f"{colors.green(model_name)}{record_repr(record, include_foreign_keys=False).lstrip(model_name)}\n"
    if record._state.db is not None and record._state.db != "default":
        msg += f"  {colors.italic('Database instance')}\n"
        msg += f"    slug: {record._state.db}\n"
    if model_name in {"Artifact", "Run"}:
        result = get_artifact_or_run_with_related(
            record,
            include_feature_link=True,
            include_fk=True,
            include_m2m=True,
            include_schema=True,
        )
        related_data = result.get("related_data", {})
        if model_name == "Artifact":
            tree = describe_artifact(record, related_data=related_data)
        else:
            tree = describe_run(record, related_data=related_data)
    elif model_name == "Collection":
        result = get_collection_with_related(record, include_fk=True)
        related_data = result.get("related_data", {})
        tree = describe_collection(record, related_data=related_data)
    elif model_name == "Transform":
        tree = describe_transform(record)
    else:
        tree = describe_header(record)
    return tree


def describe_sqlite(record):
    model_name = record.__class__.__name__
    msg = f"{colors.green(model_name)}{record_repr(record, include_foreign_keys=False).lstrip(model_name)}\n"
    if record._state.db is not None and record._state.db != "default":
        msg += f"  {colors.italic('Database instance')}\n"
        msg += f"    slug: {record._state.db}\n"

    fields = record._meta.fields
    direct_fields = []
    foreign_key_fields = []
    for f in fields:
        if f.is_relation:
            foreign_key_fields.append(f.name)
        else:
            direct_fields.append(f.name)
    if not record._state.adding:
        # prefetch foreign key relationships
        record = (
            record.__class__.objects.using(record._state.db)
            .select_related(*foreign_key_fields)
            .get(id=record.id)
        )
        # prefetch m-2-m relationships
        many_to_many_fields = []
        if model_name in {"Artifact", "Collection"}:
            many_to_many_fields.append("input_of_runs")
        if model_name == "Artifact":
            many_to_many_fields.append("feature_sets")
        record = (
            record.__class__.objects.using(record._state.db)
            .prefetch_related(*many_to_many_fields)
            .get(id=record.id)
        )
    if model_name in {"Artifact", "Run"}:
        if model_name == "Artifact":
            tree = describe_artifact(record)
        else:
            tree = describe_run(record)
    elif model_name == "Collection":
        tree = describe_collection(record)
    elif model_name == "Transform":
        tree = describe_transform(record)
    else:
        tree = describe_header(record)
    return tree


def describe_postgres_sqlite(record, return_str: bool = False) -> str | None:
    from ._describe import format_rich_tree

    if (
        not record._state.adding
        and connections[record._state.db].vendor == "postgresql"
    ):
        tree = describe_postgres(record)
    else:
        tree = describe_sqlite(record)

    return format_rich_tree(tree, return_str=return_str)
