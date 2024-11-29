import datetime
from typing import Any

from lnschema_core.models import Artifact, Collection
from rich.table import Column, Table
from rich.text import Text
from rich.tree import Tree


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


def describe_header(self: Artifact | Collection) -> Tree:
    # initialize tree
    suffix = f" - {self.suffix.removeprefix('.')}" if self.suffix else ""
    tree = Tree(Text(f"Artifact{suffix}", style="bold"), guide_style="dim")
    return tree


def describe_general(self: Artifact | Collection, tree: Tree | None = None) -> Tree:
    if tree is None:
        tree = describe_header(self)

    # add general information
    general = tree.add(Text("General", style="bold bright_cyan"))
    general.add(f".uid = '{self.uid}'")
    if hasattr(self, "n_observations") and self.n_observations:
        general.add(Text(f".n_observations = {self.n_observations}"))
    if hasattr(self, "key") and self.key:
        general.add(
            f".key (virtual) = '{self.key}'"
            if self._key_is_virtual
            else f".key = {self.key}"
        )
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
    general.add(Text.assemble(".created_at = ", highlight_time(str(self.created_at))))
    if self.transform:
        general.add(
            Text(
                f".transform = '{self.transform.name}'",
                style="khaki1",
            )
        )

    return tree


def describe_labels(
    self: Artifact | Collection,
    labels_data: dict | None = None,
    print_types: bool = False,
    tree: Tree | None = None,
    as_subtree: bool = False,
):
    from ._label_manager import (
        _get_labels,
        _get_labels_postgres,
        _print_values,
        connections,
        get_name_field,
        get_related_model,
    )

    if not self._state.adding and connections[self._state.db].vendor == "postgresql":
        labels_data = _get_labels_postgres(self, labels_data)
    if not labels_data:
        labels_data = _get_labels(self, instance=self._state.db)

    if tree is None:
        tree = describe_header(self)
    if not labels_data:
        return tree

    labels_table = Table(
        Column("Name", style="", no_wrap=True, width=NAME_WIDTH),
        Column("Type", style="dim", no_wrap=True, width=TYPE_WIDTH),
        Column("Values", width=VALUES_WIDTH, no_wrap=True),
        show_header=False,
        box=None,
        pad_edge=False,
    )
    for related_name, labels in labels_data.items():
        if not labels or related_name == "feature_sets":
            continue
        if isinstance(labels, dict):  # postgres, labels are a dict[id, name]
            print_values = _print_values(labels.values(), n=10)
        else:  # labels are a QuerySet
            field = get_name_field(labels)
            print_values = _print_values(labels.values_list(field, flat=True), n=10)
        if print_values:
            related_model = get_related_model(self, related_name)
            type_str = related_model.__get_name_with_schema__()
            labels_table.add_row(
                f".{related_name}", Text(type_str, style="dim"), print_values
            )

    if as_subtree:
        return labels_table
    else:
        labels_tree = tree.add(Text("Labels", style="bold pale_green3"))
        labels_tree.add(labels_table)
        return tree


def describe_features(
    self: Artifact | Collection,
    related_data: dict | None = None,
    print_types: bool = False,
    to_dict: bool = False,
    print_params: bool = False,
    tree: Tree | None = None,
    with_labels: bool = False,
):
    from lamindb._from_values import _print_values

    from ._feature_manager import (
        FeatureSet,
        _get_categoricals,
        _get_categoricals_postgres,
        _get_featuresets_postgres,
        _get_non_categoricals,
        connections,
        get_feature_set_by_slot_,
        get_name_field,
    )

    if tree is None:
        tree = describe_header(self)

    dataset_tree = tree.add(Text("Dataset", style="bold cornflower_blue"))
    dictionary: dict[str, Any] = {}

    if self._state.adding:
        return dictionary if to_dict else tree

    # feature sets
    feature_set_data: dict[str, tuple[str, list[str]]] = {}
    feature_data: dict[str, tuple[str, list[str]]] = {}
    if not print_params and not to_dict:
        if self.id is not None and connections[self._state.db].vendor == "postgresql":
            fs_data = _get_featuresets_postgres(self, related_data=related_data)
            for fs_id, (slot, data) in fs_data.items():
                for registry_str, feature_names in data.items():
                    feature_set = FeatureSet.get(id=fs_id)
                    feature_set_data[slot] = (feature_set, feature_names)
                    for feature_name in feature_names:
                        feature_data[feature_name] = (slot, registry_str)
        else:
            for slot, feature_set in get_feature_set_by_slot_(self).items():
                features = feature_set.members
                # features.first() is a lot slower than features[0] here
                name_field = get_name_field(features[0])
                feature_names = list(features.values_list(name_field, flat=True)[:10])
                feature_set_data[slot] = (feature_set, feature_names)
                for feature_name in feature_names:
                    feature_data[feature_name] = (slot, feature_set.registry)

    internal_feature_names: set[str] = {}  # type: ignore
    if isinstance(self, Artifact):
        feature_set = self.feature_sets.filter(registry="Feature").one_or_none()
        if feature_set is not None:
            internal_feature_names = set(
                feature_set.members.values_list("name", flat=True)
            )  # type: ignore

    # categorical feature values
    # Get the categorical data using the appropriate method
    if not self._state.adding and connections[self._state.db].vendor == "postgresql":
        categoricals = _get_categoricals_postgres(
            self,
            related_data=related_data,
            print_params=print_params,
        )
    else:
        categoricals = _get_categoricals(
            self,
            print_params=print_params,
        )

    # Get non-categorical features
    non_categoricals = _get_non_categoricals(
        self,
        print_params=print_params,
    )

    # Combine all features
    internal_data = {}
    external_data = []

    # Process all Features and sort into internal/external
    for features, is_list_type in [(categoricals, False), (non_categoricals, True)]:
        for (feature_name, feature_dtype), values in sorted(features.items()):
            # Handle dictionary conversion - no separation needed for dictionary
            if to_dict:
                dict_value = values if len(values) > 1 else next(iter(values))
                dictionary[feature_name] = dict_value
            else:
                # Format message
                printed_values = (
                    _print_values(sorted(values), n=10, quotes=False)
                    if not is_list_type or not feature_dtype.startswith("list")
                    else sorted(values)
                )

                # Sort into internal/external
                if feature_name in internal_feature_names:
                    internal_data[feature_name] = (
                        feature_name,
                        Text(feature_dtype, style="dim"),
                        printed_values,
                    )
                else:
                    external_data.append(
                        (
                            feature_name,
                            Text(feature_dtype, style="dim"),
                            printed_values,
                        )
                    )

    if to_dict:
        return dictionary

    # Internal features from the Feature registry
    internal_features_slot: dict[str, list] = {}
    for feature_name, feature_row in internal_data.items():
        slot, registry_str = feature_data.get(feature_name)
        if slot not in internal_features_slot:
            internal_features_slot[slot] = []
        else:
            internal_features_slot[slot].append(feature_row)

    for slot, feature_rows in internal_features_slot.items():
        feature_set = feature_set_data[slot][0]
        dataset_tree.add(
            _create_feature_table(
                Text.assemble(
                    (slot, "violet"),
                    (f" ← {feature_set.n} [{feature_set.registry}]", "dim"),
                ),
                feature_rows,
            )
        )

    # Internal features from other registries (e.g. bionty.Gene)
    for slot, (feature_set, feature_names) in feature_set_data.items():
        if feature_set.registry == "Feature":
            continue
        features = [
            (feature_name, Text(str(feature_set.dtype), style="dim"), "")
            for feature_name in feature_names
        ]
        dataset_tree.add(
            _create_feature_table(
                Text.assemble(
                    (slot, "violet"),
                    (f" ← {feature_set.n} [{feature_set.registry}]", "dim"),
                ),
                features,
            )
        )

    # Annotations (external features)
    annotations_tree = tree.add(Text("Annotations", style="bold dark_orange"))
    annotations_tree.add(
        _create_feature_table(Text.assemble(("Features", "pale_green3")), external_data)
    )

    if with_labels:
        labels_table = describe_labels(self, as_subtree=True)
        annotations_tree.add(Text("Labels", style="bold pale_green3"))
        annotations_tree.add(labels_table)
    return tree


# Create combined tables for each feature group that include the header
def _create_feature_table(name: str, data: list) -> Table:
    table = Table(
        Column(name, style="", no_wrap=True, width=NAME_WIDTH),
        Column("", style="dim", no_wrap=True, width=TYPE_WIDTH),
        Column("", width=VALUES_WIDTH, no_wrap=True),
        show_header=True,
        box=None,
        pad_edge=False,
    )
    for row in data:
        table.add_row(*row)
    return table
