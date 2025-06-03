# ruff: noqa: TC004
from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from datetime import date, datetime
from itertools import compress
from typing import TYPE_CHECKING, Any, MutableMapping

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData
from django.contrib.postgres.aggregates import ArrayAgg
from django.db import connections
from django.db.models import Aggregate, ProtectedError, Subquery
from lamin_utils import logger
from lamindb_setup.core.hashing import hash_set
from lamindb_setup.core.upath import create_path
from rich.table import Column, Table
from rich.text import Text

from lamindb.core.storage import LocalPathClasses
from lamindb.errors import DoesNotExist, ValidationError
from lamindb.models._from_values import _format_values
from lamindb.models.feature import (
    serialize_pandas_dtype,
    suggest_categorical_for_str_iterable,
)
from lamindb.models.save import save
from lamindb.models.schema import DICT_KEYS_TYPE, Schema
from lamindb.models.sqlrecord import (
    REGISTRY_UNIQUE_FIELD,
    get_name_field,
    transfer_fk_to_default_db_bulk,
    transfer_to_default_db,
)

from ..base import deprecated
from ._describe import (
    NAME_WIDTH,
    TYPE_WIDTH,
    VALUES_WIDTH,
    describe_header,
    format_rich_tree,
)
from ._django import get_artifact_with_related
from ._label_manager import _get_labels, describe_labels
from ._relations import (
    dict_related_model_to_related_name,
)
from .feature import Feature, FeatureValue, parse_dtype
from .run import FeatureManager, FeatureManagerRun, Run
from .sqlrecord import SQLRecord
from .ulabel import ULabel

if TYPE_CHECKING:
    from rich.tree import Tree

    from lamindb.base.types import FieldAttr
    from lamindb.models import (
        Artifact,
        Collection,
        IsLink,
    )
    from lamindb.models.query_set import QuerySet


class FeatureManagerArtifact(FeatureManager):
    """Feature manager."""

    pass


def get_accessor_by_registry_(host: Artifact | Collection) -> dict:
    dictionary = {
        field.related_model.__get_name_with_module__(): field.name
        for field in host._meta.related_objects
    }
    dictionary["Feature"] = "features"
    dictionary["ULabel"] = "ulabels"
    return dictionary


def get_schema_by_slot_(host: Artifact) -> dict:
    # if the host is not yet saved
    if host._state.adding:
        if hasattr(host, "_staged_feature_sets"):
            return host._staged_feature_sets
        else:
            return {}
    host_db = host._state.db
    kwargs = {"artifact_id": host.id}
    # otherwise, we need a query
    links_schema = (
        host.feature_sets.through.objects.using(host_db)
        .filter(**kwargs)
        .select_related("schema")
    )
    return {fsl.slot: fsl.schema for fsl in links_schema}


def get_label_links(
    host: Artifact | Collection, registry: str, feature: Feature
) -> QuerySet:
    kwargs = {"artifact_id": host.id, "feature_id": feature.id}
    link_records = (
        getattr(host, host.features._accessor_by_registry[registry])  # type: ignore
        .through.objects.using(host._state.db)
        .filter(**kwargs)
    )
    return link_records


def get_schema_links(host: Artifact | Collection) -> QuerySet:
    kwargs = {"artifact_id": host.id}
    links_schema = host.feature_sets.through.objects.filter(**kwargs)
    return links_schema


def get_link_attr(link: IsLink | type[IsLink], data: Artifact | Collection) -> str:
    link_model_name = link.__class__.__name__
    if link_model_name in {"Registry", "ModelBase"}:  # we passed the type of the link
        link_model_name = link.__name__  # type: ignore
    return link_model_name.replace(data.__class__.__name__, "").lower()


# Custom aggregation for SQLite
class GroupConcat(Aggregate):
    function = "GROUP_CONCAT"
    template = '%(function)s(%(expressions)s, ", ")'


def custom_aggregate(field, using: str):
    if connections[using].vendor == "postgresql":
        return ArrayAgg(field)
    else:
        return GroupConcat(field)


def _get_categoricals_postgres(
    self: Artifact | Collection | Run,
    related_data: dict | None = None,
) -> dict[tuple[str, str], set[str]]:
    """Get categorical features and their values using PostgreSQL-specific optimizations."""
    if not related_data:
        if self.__class__.__name__ == "Artifact":
            artifact_meta = get_artifact_with_related(
                self, include_feature_link=True, include_m2m=True
            )
            related_data = artifact_meta.get("related_data", {})
        else:
            related_data = {}

    # Process m2m data
    m2m_data = related_data.get("m2m", {}) if related_data else {}
    m2m_name = {}
    for related_name, values in m2m_data.items():
        link_model = getattr(self.__class__, related_name).through
        related_model_name = link_model.__name__.replace(
            self.__class__.__name__, ""
        ).lower()
        m2m_name[related_model_name] = values

    # Get feature information
    links_data = related_data.get("link", {}) if related_data else {}
    feature_dict = {
        id: (name, dtype)
        for id, name, dtype in Feature.objects.using(self._state.db).values_list(
            "id", "name", "dtype"
        )
    }

    # Build result dictionary
    result = defaultdict(set)
    for link_name, link_values in links_data.items():
        related_name = link_name.removeprefix("links_").replace("_", "")
        if not link_values:
            continue

        for link_value in link_values:
            feature_id = link_value.get("feature")
            if feature_id is None:
                continue

            feature_name, feature_dtype = feature_dict.get(feature_id)
            label_id = link_value.get(related_name)
            label_name = m2m_name.get(related_name, {}).get(label_id)
            if label_name:
                result[(feature_name, feature_dtype)].add(label_name)

    return dict(result)


def _get_categoricals(
    self: Artifact | Collection,
) -> dict[tuple[str, str], set[str]]:
    """Get categorical features and their values using the default approach."""
    result = defaultdict(set)
    for _, links in _get_labels(self, links=True, instance=self._state.db).items():
        for link in links:
            if hasattr(link, "feature_id") and link.feature_id is not None:
                feature = Feature.objects.using(self._state.db).get(id=link.feature_id)
                link_attr = get_link_attr(link, self)
                label = getattr(link, link_attr)
                name_attr = (
                    "name" if hasattr(label, "name") else label.__class__._name_field
                )
                label_name = getattr(label, name_attr)
                result[(feature.name, feature.dtype)].add(label_name)

    return dict(result)


def _get_non_categoricals(
    self,
) -> dict[tuple[str, str], set[Any]]:
    """Get non-categorical features and their values."""
    from .artifact import Artifact
    from .run import Run

    non_categoricals = {}

    if self.id is not None and isinstance(self, (Artifact, Run)):
        attr_name = "feature"
        _feature_values = (
            getattr(self, f"_{attr_name}_values")
            .values(f"{attr_name}__name", f"{attr_name}__dtype")
            .annotate(values=custom_aggregate("value", self._state.db))
            .order_by(f"{attr_name}__name")
        )

        for fv in _feature_values:
            feature_name = fv[f"{attr_name}__name"]
            feature_dtype = fv[f"{attr_name}__dtype"]
            values = fv["values"]

            # Convert single values to sets
            if not isinstance(values, (list, dict, set)):
                values = {values}
            elif (
                isinstance(values, list)
                and feature_dtype != "dict"
                and not feature_dtype.startswith("list")
            ):
                try:
                    values = set(values)
                except TypeError:
                    # TypeError: unhashable type: 'list' if values is list[list]
                    pass

            # Handle special datetime types
            if feature_dtype == "datetime":
                values = {datetime.fromisoformat(value) for value in values}
            if feature_dtype == "date":
                values = {date.fromisoformat(value) for value in values}

            non_categoricals[(feature_name, feature_dtype)] = values

    return non_categoricals


def _get_schemas_postgres(
    self: Artifact | Collection,
    related_data: dict | None = None,
) -> dict:
    if not related_data:
        artifact_meta = get_artifact_with_related(self, include_schema=True)
        related_data = artifact_meta.get("related_data", {})

    fs_data = related_data.get("schemas", {}) if related_data else {}
    return fs_data


def _create_feature_table(
    name: str, registry_str: str, data: list, show_header: bool = False
) -> Table:
    """Create a Rich table for a feature group."""
    table = Table(
        Column(name, style="", no_wrap=True, width=NAME_WIDTH),
        Column(registry_str, style="dim", no_wrap=True, width=TYPE_WIDTH),
        Column("", width=VALUES_WIDTH, no_wrap=True),
        show_header=show_header,
        box=None,
        pad_edge=False,
    )
    for row in data:
        table.add_row(*row)
    return table


def describe_features(
    self: Artifact,
    related_data: dict | None = None,
    to_dict: bool = False,
    tree: Tree | None = None,
    with_labels: bool = False,
):
    """Describe features of an artifact or collection."""
    from .artifact import Artifact

    # initialize tree
    if tree is None:
        tree = describe_header(self)

    dictionary: dict[str, Any] = {}

    if self._state.adding:
        return dictionary if to_dict else tree

    # feature sets
    schema_data: dict[str, tuple[str, list[str]]] = {}
    feature_data: dict[str, tuple[str, list[str]]] = {}
    if not to_dict:
        if self.id is not None and connections[self._state.db].vendor == "postgresql":
            fs_data = _get_schemas_postgres(self, related_data=related_data)
            for fs_id, (slot, data) in fs_data.items():
                for registry_str, feature_names in data.items():
                    # prevent projects show up as features
                    if registry_str == "Project":
                        continue
                    schema = Schema.objects.using(self._state.db).get(id=fs_id)
                    schema_data[slot] = (schema, feature_names)
                    for feature_name in feature_names:
                        feature_data[feature_name] = (slot, registry_str)
            schema_data.update(
                {
                    slot: (schema, schema.n)
                    for slot, schema in get_schema_by_slot_(self).items()
                    if slot not in schema_data
                }
            )
        else:
            for slot, schema in get_schema_by_slot_(self).items():
                features = schema.members
                if features.exists():
                    # features.first() is a lot slower than features[0] here
                    name_field = get_name_field(features[0])
                    feature_names = list(
                        features.values_list(name_field, flat=True)[:20]
                    )
                    schema_data[slot] = (schema, feature_names)
                    for feature_name in feature_names:
                        feature_data[feature_name] = (slot, schema.itype)
                else:
                    schema_data[slot] = (schema, schema.n)

    internal_feature_names: dict[str, str] = {}
    if isinstance(self, Artifact):
        feature_sets = self.feature_sets.filter(itype="Feature").all()
        internal_feature_names = {}
        if len(feature_sets) > 0:
            for schema in feature_sets:
                internal_feature_names.update(
                    dict(schema.members.values_list("name", "dtype"))
                )

    # categorical feature values
    # Get the categorical data using the appropriate method
    if not self._state.adding and connections[self._state.db].vendor == "postgresql":
        categoricals = _get_categoricals_postgres(
            self,
            related_data=related_data,
        )
    else:
        categoricals = _get_categoricals(
            self,
        )

    # Get non-categorical features
    non_categoricals = _get_non_categoricals(
        self,
    )

    # Process all Features containing labels and sort into internal/external
    internal_feature_labels = {}
    external_data = []
    for features, is_list_type in [(categoricals, False), (non_categoricals, True)]:
        for (feature_name, feature_dtype), values in sorted(features.items()):
            # Handle dictionary conversion
            if to_dict:
                dict_value = values if len(values) > 1 else next(iter(values))
                dictionary[feature_name] = dict_value
                continue

            # Format message
            printed_values = (
                _format_values(sorted(values), n=10, quotes=False)
                if not is_list_type or not feature_dtype.startswith("list")
                else sorted(values)
            )

            # Sort into internal/external
            feature_info = (
                feature_name,
                Text(feature_dtype, style="dim"),
                printed_values,
            )
            if feature_name in internal_feature_names:
                internal_feature_labels[feature_name] = feature_info
            else:
                external_data.append(feature_info)

    if to_dict:
        return dictionary

    # Dataset features section
    # internal features that contain labels (only `Feature` features contain labels)
    internal_feature_labels_slot: dict[str, list] = {}
    for feature_name, feature_row in internal_feature_labels.items():
        slot, _ = feature_data.get(feature_name)
        internal_feature_labels_slot.setdefault(slot, []).append(feature_row)

    int_features_tree_children = []
    for slot, (schema, feature_names_or_n) in schema_data.items():
        if isinstance(feature_names_or_n, int):
            feature_rows = []
        else:
            feature_names = feature_names_or_n
            if slot in internal_feature_labels_slot:
                # add internal Feature features with labels
                feature_rows = internal_feature_labels_slot[slot]
                # add internal Feature features without labels
                feature_rows += [
                    (
                        feature_name,
                        Text(
                            str(internal_feature_names.get(feature_name)), style="dim"
                        ),
                        "",
                    )
                    for feature_name in feature_names
                    if feature_name and feature_name not in internal_feature_labels
                ]
            else:
                # add internal non-Feature features without labels
                feature_rows = [
                    (
                        feature_name,
                        Text(
                            str(
                                internal_feature_names.get(feature_name)
                                if feature_name in internal_feature_names
                                else schema.dtype
                            ),
                            style="dim",
                        ),
                        "",
                    )
                    for feature_name in feature_names
                    if feature_name
                ]
        int_features_tree_children.append(
            _create_feature_table(
                Text.assemble(
                    (slot, "violet"),
                    (" • ", "dim"),
                    (str(schema.n), "pink1"),
                ),
                Text.assemble((f"[{schema.itype}]", "pink1")),
                feature_rows,
                show_header=True,
            )
        )
    ## internal features from the non-`Feature` registry
    if int_features_tree_children:
        dataset_tree = tree.add(
            Text.assemble(
                ("Dataset features", "bold bright_magenta"),
            )
        )
        for child in int_features_tree_children:
            dataset_tree.add(child)

    # Linked features
    ext_features_tree_children = []
    if external_data:
        ext_features_tree_children.append(
            _create_feature_table(
                "",
                "",
                external_data,
            )
        )
    # ext_features_tree = None
    ext_features_header = Text("Linked features", style="bold dark_orange")
    if ext_features_tree_children:
        ext_features_tree = tree.add(ext_features_header)
        for child in ext_features_tree_children:
            ext_features_tree.add(child)
    if with_labels:
        # avoid querying the db if the labels were queried already
        labels_data = related_data.get("m2m") if related_data is not None else None
        labels_tree = describe_labels(self, labels_data=labels_data, as_subtree=True)
        if labels_tree:
            tree.add(labels_tree)

    return tree


def is_valid_datetime_str(date_string: str) -> bool | str:
    try:
        dt = datetime.fromisoformat(date_string)
        return dt.isoformat()
    except ValueError:
        return False


def infer_feature_type_convert_json(
    key: str, value: Any, mute: bool = False, str_as_ulabel: bool = True
) -> tuple[str, Any, str]:
    message = ""
    if isinstance(value, bool):
        return "bool", value, message
    elif isinstance(value, int):
        return "int", value, message
    elif isinstance(value, float):
        return "float", value, message
    elif isinstance(value, date):
        return "date", value.isoformat(), message
    elif isinstance(value, datetime):
        return "datetime", value.isoformat(), message
    elif isinstance(value, str):
        if datetime_str := is_valid_datetime_str(value):
            dt_type = (
                "date" if len(value) == 10 else "datetime"
            )  # YYYY-MM-DD is exactly 10 characters
            sanitized_value = datetime_str[:10] if dt_type == "date" else datetime_str  # type: ignore
            return dt_type, sanitized_value, message  # type: ignore
        else:
            return "cat ? str", value, message
    elif isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        if isinstance(value, (pd.Series, np.ndarray, pd.Categorical)):
            dtype = serialize_pandas_dtype(value.dtype)
            if dtype == "str":
                # ndarray doesn't know categorical, so there was no conscious choice
                # offer both options
                if isinstance(value, np.ndarray):
                    dtype = "cat ? str"
                else:
                    # suggest to create a categorical if there are few unique values
                    message = suggest_categorical_for_str_iterable(value, key)
                    if message:
                        message = f"  # {message}"
            return dtype, list(value), message
        if isinstance(value, dict):
            return "dict", value, message
        if len(value) > 0:  # type: ignore
            first_element_type = type(next(iter(value)))
            if all(isinstance(elem, first_element_type) for elem in value):
                if first_element_type is bool:
                    return "list[bool]", value, message
                elif first_element_type is int:
                    return "list[int]", value, message
                elif first_element_type is float:
                    return "list[float]", value, message
                elif first_element_type is str:
                    return ("list[cat ? str]", value, message)
                elif first_element_type == SQLRecord:
                    return (
                        f"list[cat[{first_element_type.__get_name_with_module__()}]]",
                        value,
                        message,
                    )
    elif isinstance(value, SQLRecord):
        return (f"cat[{value.__class__.__get_name_with_module__()}]", value, message)
    if not mute:
        logger.warning(f"cannot infer feature type of: {value}, returning '?")
    return "?", value, message


def __init__(self, host: Artifact | Collection | Run):
    self._host = host
    self._slots = None
    self._accessor_by_registry_ = None


def __repr__(self) -> str:
    return describe(self, return_str=True)  # type: ignore


def describe(self, return_str: bool = False) -> str | None:
    tree = describe_features(self._host)  # type: ignore
    return format_rich_tree(tree, fallback="no linked features", return_str=return_str)


def get_values(self) -> dict[str, Any]:
    """Get feature values as a dictionary."""
    return describe_features(self._host, to_dict=True)  # type: ignore


@deprecated("slots[slot].members")
def __getitem__(self, slot) -> QuerySet:
    if slot not in self.slots:
        raise ValueError(
            f"No linked feature set for slot: {slot}\nDid you get validation"
            " warnings? Only features that match registered features get validated"
            " and linked."
        )
    schema = self.slots[slot]
    orm_name = schema.itype
    return getattr(schema, self._accessor_by_registry[orm_name]).all()


def filter_base(cls, _skip_validation: bool = True, **expression) -> QuerySet:
    from .artifact import Artifact

    model = Feature
    value_model = FeatureValue
    keys_normalized = [key.split("__")[0] for key in expression]
    if not _skip_validation:
        validated = model.validate(keys_normalized, field="name", mute=True)
        if sum(validated) != len(keys_normalized):
            raise ValidationError(
                f"Some keys in the filter expression are not registered as features: {np.array(keys_normalized)[~validated]}"
            )
    new_expression = {}
    features = model.filter(name__in=keys_normalized).all().distinct()
    feature_param = "feature"
    for key, value in expression.items():
        split_key = key.split("__")
        normalized_key = split_key[0]
        comparator = ""
        if len(split_key) == 2:
            comparator = f"__{split_key[1]}"
        feature = features.get(name=normalized_key)
        # non-categorical features
        if not feature.dtype.startswith("cat") and not feature.dtype.startswith(
            "list[cat"
        ):
            if comparator == "__isnull":
                if cls == FeatureManagerArtifact:
                    from .artifact import ArtifactFeatureValue

                    if value:  # True
                        return Artifact.objects.exclude(
                            id__in=Subquery(
                                ArtifactFeatureValue.objects.filter(
                                    featurevalue__feature=feature
                                ).values("artifact_id")
                            )
                        )
                    else:
                        return Artifact.objects.exclude(
                            id__in=Subquery(
                                ArtifactFeatureValue.objects.filter(
                                    featurevalue__feature=feature
                                ).values("artifact_id")
                            )
                        )
            if comparator in {"__startswith", "__contains"}:
                logger.important(
                    f"currently not supporting `{comparator}`, using `__icontains` instead"
                )
                comparator = "__icontains"
            expression = {feature_param: feature, f"value{comparator}": value}
            feature_values = value_model.filter(**expression)
            new_expression[f"_{feature_param}_values__id__in"] = feature_values
        # categorical features
        elif isinstance(value, (str, SQLRecord, bool)):
            if comparator == "__isnull":
                if cls == FeatureManagerArtifact:
                    result = parse_dtype(feature.dtype)[0]
                    kwargs = {
                        f"links_{result['registry'].__name__.lower()}__feature": feature
                    }
                    if value:  # True
                        return Artifact.objects.exclude(**kwargs)
                    else:
                        return Artifact.objects.filter(**kwargs)
            else:
                # because SQL is sensitive to whether querying with __in or not
                # and might return multiple equivalent records for the latter
                # we distinguish cases in which we have multiple label matches vs. one
                label = None
                labels = None
                result = parse_dtype(feature.dtype)[0]
                label_registry = result["registry"]
                if isinstance(value, str):
                    field_name = result["field"].field.name
                    # we need the comparator here because users might query like so
                    # ln.Artifact.filter(experiment__contains="Experi")
                    expression = {f"{field_name}{comparator}": value}
                    labels = result["registry"].filter(**expression).all()
                    if len(labels) == 0:
                        raise DoesNotExist(
                            f"Did not find a {label_registry.__name__} matching `{field_name}{comparator}={value}`"
                        )
                    elif len(labels) == 1:
                        label = labels[0]
                elif isinstance(value, SQLRecord):
                    label = value
                accessor_name = (
                    label_registry.artifacts.through.artifact.field._related_name
                )
                new_expression[f"{accessor_name}__feature"] = feature
                if label is not None:
                    # simplified query if we have exactly one label
                    new_expression[
                        f"{accessor_name}__{label_registry.__name__.lower()}"
                    ] = label
                else:
                    new_expression[
                        f"{accessor_name}__{label_registry.__name__.lower()}__in"
                    ] = labels
            # if passing a list of records, we want to
            # find artifacts that are annotated by all of them at the same
            # time; hence, we don't want the __in construct that we use to match strings
            # https://laminlabs.slack.com/archives/C04FPE8V01W/p1688328084810609
    if not (new_expression):
        raise NotImplementedError
    if cls == FeatureManagerArtifact:
        return Artifact.objects.filter(**new_expression)
    elif cls == FeatureManagerRun:
        return Run.objects.filter(**new_expression)


@classmethod  # type: ignore
@deprecated("the filter() registry classmethod")
def filter(cls, **expression) -> QuerySet:
    """Query artifacts by features."""
    return filter_base(cls, _skip_validation=False, **expression)


@classmethod  # type: ignore
@deprecated("the filter() registry classmethod")
def get(cls, **expression) -> SQLRecord:
    """Query a single artifact by feature."""
    return filter_base(cls, _skip_validation=False, **expression).one()


@property  # type: ignore
def slots(self) -> dict[str, Schema]:
    """Schema by slot.

    Example::

        artifact.features.slots
        #> {'var': <Schema: var>, 'obs': <Schema: obs>}
    """
    if self._slots is None:
        self._slots = get_schema_by_slot_(self._host)
    return self._slots


@property  # type: ignore
def _accessor_by_registry(self):
    """Accessor by ORM."""
    if self._accessor_by_registry_ is None:
        self._accessor_by_registry_ = get_accessor_by_registry_(self._host)
    return self._accessor_by_registry_


def add_label_feature_links(
    self,
    features_labels,
    *,
    label_ref_is_name: bool | None = None,
    feature_ref_is_name: bool | None = None,
):
    if list(features_labels.keys()) != ["ULabel"]:
        related_names = dict_related_model_to_related_name(self._host.__class__)
    else:
        related_names = {"ULabel": "ulabels"}
    for class_name, registry_features_labels in features_labels.items():
        related_name = related_names[class_name]  # e.g., "ulabels"
        IsLink = getattr(self._host, related_name).through
        field_name = f"{get_link_attr(IsLink, self._host)}_id"  # e.g., ulabel_id
        links = [
            IsLink(
                **{
                    "artifact_id": self._host.id,
                    "feature_id": feature.id,
                    field_name: label.id,
                    "feature_ref_is_name": feature_ref_is_name,
                    "label_ref_is_name": label_ref_is_name,
                }
            )
            for (feature, label) in registry_features_labels
        ]
        # a link might already exist
        try:
            save(links, ignore_conflicts=False)
        except Exception:
            save(links, ignore_conflicts=True)
        # now delete links that were previously saved without a feature
        IsLink.filter(
            **{
                "artifact_id": self._host.id,
                "feature_id": None,
                f"{field_name}__in": [l.id for _, l in registry_features_labels],
            }
        ).all().delete()


def _add_values(
    self,
    values: dict[str, str | int | float | bool],
    feature_param_field: FieldAttr,
    str_as_ulabel: bool = True,
) -> None:
    """Curate artifact with features & values.

    Args:
        values: A dictionary of keys (features) & values (labels, numbers, booleans).
        feature_param_field: The field of a reference registry to map keys of the
            dictionary.
    """
    from .._tracked import get_current_tracked_run

    # rename to distinguish from the values inside the dict
    dictionary = values
    keys = dictionary.keys()
    if isinstance(keys, DICT_KEYS_TYPE):
        keys = list(keys)  # type: ignore
    # deal with other cases later
    assert all(isinstance(key, str) for key in keys)  # noqa: S101
    registry = feature_param_field.field.model
    value_model = FeatureValue
    model_name = "Feature"
    records = registry.from_values(keys, field=feature_param_field, mute=True)
    if len(records) != len(keys):
        not_validated_keys = [key for key in keys if key not in records.list("name")]
        not_validated_keys_dtype_message = [
            (key, infer_feature_type_convert_json(key, dictionary[key]))
            for key in not_validated_keys
        ]
        run = get_current_tracked_run()
        if run is not None:
            name = f"{run.transform.type}[{run.transform.key}]"
            type_hint = f"""  {model_name.lower()}_type = ln.{model_name}(name='{name}', is_type=True).save()"""
            elements = [type_hint]
            type_kwarg = f", type={model_name.lower()}_type"
        else:
            elements = []
            type_kwarg = ""
        elements += [
            f"  ln.{model_name}(name='{key}', dtype='{dtype}'{type_kwarg}).save(){message}"
            for key, (dtype, _, message) in not_validated_keys_dtype_message
        ]
        hint = "\n".join(elements)
        msg = (
            f"These keys could not be validated: {not_validated_keys}\n"
            f"Here is how to create a {model_name.lower()}:\n\n{hint}"
        )
        raise ValidationError(msg)

    # figure out which of the values go where
    features_labels = defaultdict(list)
    _feature_values = []
    not_validated_values = []
    for feature in records:
        value = dictionary[feature.name]
        inferred_type, converted_value, _ = infer_feature_type_convert_json(
            feature.name,
            value,
            mute=True,
            str_as_ulabel=str_as_ulabel,
        )
        if feature.dtype == "num":
            if inferred_type not in {"int", "float"}:
                raise TypeError(
                    f"Value for feature '{feature.name}' with type {feature.dtype} must be a number"
                )
        elif feature.dtype.startswith("cat"):
            if inferred_type != "?":
                if not (
                    inferred_type.startswith("cat") or isinstance(value, SQLRecord)
                ):
                    raise TypeError(
                        f"Value for feature '{feature.name}' with type '{feature.dtype}' must be a string or record."
                    )
        elif (feature.dtype == "str" and feature.dtype not in inferred_type) or (
            feature.dtype != "str" and feature.dtype != inferred_type
        ):
            raise ValidationError(
                f"Expected dtype for '{feature.name}' is '{feature.dtype}', got '{inferred_type}'"
            )
        if not feature.dtype.startswith("cat"):
            filter_kwargs = {model_name.lower(): feature, "value": converted_value}
            feature_value, _ = value_model.get_or_create(**filter_kwargs)
            _feature_values.append(feature_value)
        else:
            if isinstance(value, SQLRecord) or (
                isinstance(value, Iterable) and isinstance(next(iter(value)), SQLRecord)
            ):
                if isinstance(value, SQLRecord):
                    label_records = [value]
                else:
                    label_records = value  # type: ignore
                for record in label_records:
                    if record._state.adding:
                        raise ValidationError(
                            f"Please save {record} before annotation."
                        )
                    features_labels[record.__class__.__get_name_with_module__()].append(
                        (feature, record)
                    )
            else:
                if isinstance(value, str):
                    values = [value]  # type: ignore
                else:
                    values = value  # type: ignore
                if "ULabel" not in feature.dtype:
                    feature.dtype += "[ULabel]"
                    feature.save()
                validated = ULabel.validate(values, field=ULabel.name, mute=True)
                values_array = np.array(values)
                validated_values = values_array[validated]
                if validated.sum() != len(values):
                    not_validated_values += values_array[~validated].tolist()
                label_records = ULabel.from_values(
                    validated_values, field=ULabel.name, mute=True
                )  # type: ignore
                features_labels["ULabel"] += [
                    (feature, label_record) for label_record in label_records
                ]
    if not_validated_values:
        not_validated_values.sort()
        hint = f"  ulabels = ln.ULabel.from_values({not_validated_values}, create=True).save()\n"
        msg = (
            f"These values could not be validated: {not_validated_values}\n"
            f"Here is how to create ulabels for them:\n\n{hint}"
        )
        raise ValidationError(msg)
    # TODO: create an explicit version of this
    # if not is_param:
    #     # check if _expect_many is false for _all_ records
    #     if any(record._expect_many for record in records):
    #         updated_features = []
    #         for record in records:
    #             if record._expect_many:
    #                 record._expect_many = False
    #                 record.save()
    #                 updated_features.append(record.name)
    #         if any(updated_features):
    #             logger.important(
    #                 f"changed observational unit to Artifact for features: {', '.join(updated_features)}"
    #             )
    # bulk add all links
    if features_labels:
        add_label_feature_links(self, features_labels)
    if _feature_values:
        to_insert_feature_values = [
            record for record in _feature_values if record._state.adding
        ]
        if to_insert_feature_values:
            save(to_insert_feature_values)
        dict_typed_features = [
            getattr(record, model_name.lower())
            for record in _feature_values
            if getattr(record, model_name.lower()).dtype == "dict"
        ]
        IsLink = self._host._feature_values.through
        valuefield_id = "featurevalue_id"
        host_class_lower = self._host.__class__.__get_name_with_module__().lower()
        if dict_typed_features:
            # delete all previously existing anotations with dictionaries
            kwargs = {
                f"links_{host_class_lower}__{host_class_lower}_id": self._host.id,
                f"{model_name.lower()}__in": dict_typed_features,
            }
            try:
                value_model.filter(**kwargs).all().delete()
            except ProtectedError:
                pass
        # add new feature links
        links = [
            IsLink(
                **{
                    f"{host_class_lower}_id": self._host.id,
                    valuefield_id: feature_value.id,
                }
            )
            for feature_value in _feature_values
        ]
        # a link might already exist, to avoid raising a unique constraint
        # error, ignore_conflicts
        save(links, ignore_conflicts=True)


def add_values_features(
    self,
    values: dict[str, str | int | float | bool],
    feature_field: FieldAttr = Feature.name,
    str_as_ulabel: bool = True,
) -> None:
    """Curate artifact with features & values.

    Args:
        values: A dictionary of keys (features) & values (labels, numbers, booleans).
        feature_field: The field of a reference registry to map keys of the
            dictionary.
        str_as_ulabel: Whether to interpret string values as ulabels.
    """
    _add_values(self, values, feature_field, str_as_ulabel=str_as_ulabel)


def remove_values(
    self,
    feature: str | Feature,
    *,
    value: Any | None = None,
):
    """Remove value annotations for a given feature.

    Args:
        feature: The feature for which to remove values.
        value: An optional value to restrict removal to a single value.

    """
    from .artifact import Artifact

    if isinstance(feature, str):
        feature = Feature.get(name=feature)
    filter_kwargs = {"feature": feature}
    if feature.dtype.startswith("cat["):  # type: ignore
        feature_registry = feature.dtype.replace("cat[", "").replace("]", "")  # type: ignore
        if value is not None:
            assert isinstance(value, SQLRecord)  # noqa: S101
            # the below uses our convention for field names in link models
            link_name = (
                feature_registry.split(".")[1]
                if "." in feature_registry
                else feature_registry
            ).lower()
            filter_kwargs[link_name] = value
        if feature_registry == "ULabel":
            link_attribute = "links_ulabel"
        else:
            link_models_on_models = {
                getattr(
                    Artifact, obj.related_name
                ).through.__get_name_with_module__(): obj.related_model.__get_name_with_module__()
                for obj in Artifact._meta.related_objects
                if obj.related_model.__get_name_with_module__() == feature_registry
            }
            link_attribute = {
                obj.related_name
                for obj in Artifact._meta.related_objects
                if obj.related_model.__get_name_with_module__() in link_models_on_models
            }.pop()
        getattr(self._host, link_attribute).filter(**filter_kwargs).all().delete()
    else:
        if value is not None:
            filter_kwargs["value"] = value
        feature_values = self._host._feature_values.filter(**filter_kwargs)
        self._host._feature_values.remove(*feature_values)
        # this might leave a dangling feature_value record
        # but we don't want to pay the price of making another query just to remove this annotation
        # we can clean the FeatureValue registry periodically if we want to


def _add_schema(self, schema: Schema, slot: str) -> None:
    """Annotate artifact with a schema.

    Args:
        schema: `Schema` A schema record.
        slot: `str` The slot that marks where the schema is stored in
            the artifact.
    """
    # TODO: deprecate as soon as we have the Schema-based curators
    if self._host._state.adding:
        raise ValueError(
            "Please save the artifact or collection before adding a feature set!"
        )
    host_db = self._host._state.db
    schema.save(using=host_db)
    kwargs = {
        "artifact_id": self._host.id,
        "schema": schema,
        "slot": slot,
    }
    link_record = (
        self._host.feature_sets.through.objects.using(host_db)
        .filter(**kwargs)
        .one_or_none()
    )
    if link_record is None:
        self._host.feature_sets.through(**kwargs).save(using=host_db)
        if slot in self.slots:
            logger.debug(f"replaced existing {slot} feature set")
        self._slots[slot] = schema  # type: ignore


def _unify_staged_feature_sets_by_hash(
    feature_sets: MutableMapping[str, Schema],
):
    unique_values: dict[str, Any] = {}

    for key, value in feature_sets.items():
        value_hash = value.hash  # Assuming each value has a .hash attribute
        if value_hash in unique_values:
            feature_sets[key] = unique_values[value_hash]
        else:
            unique_values[value_hash] = value

    return feature_sets


def _add_from(self, data: Artifact | Collection, transfer_logs: dict = None):
    """Transfer features from a artifact or collection."""
    # This only covers feature sets
    if transfer_logs is None:
        transfer_logs = {"mapped": [], "transferred": [], "run": None}
    from lamindb import settings

    using_key = settings._using_key
    for slot, schema in data.features.slots.items():  # type: ignore
        members = schema.members
        if len(members) == 0:
            continue
        registry = members[0].__class__
        # note here the features are transferred based on an unique field
        field = REGISTRY_UNIQUE_FIELD.get(registry.__name__.lower(), "uid")
        # this will be e.g. be a list of ontology_ids or uids
        member_uids = list(members.values_list(field, flat=True))
        validated = registry.validate(member_uids, field=field, mute=True)
        new_members_uids = list(compress(member_uids, ~validated))
        new_members = members.filter(**{f"{field}__in": new_members_uids}).all()
        n_new_members = len(new_members)
        if n_new_members > 0:
            # transfer foreign keys needs to be run before transfer to default db
            transfer_fk_to_default_db_bulk(
                new_members, using_key, transfer_logs=transfer_logs
            )
            for feature in new_members:
                # not calling save=True here as in labels, because want to
                # bulk save below
                # transfer_fk is set to False because they are already transferred
                # in the previous step transfer_fk_to_default_db_bulk
                transfer_to_default_db(
                    feature, using_key, transfer_fk=False, transfer_logs=transfer_logs
                )
            logger.info(f"saving {n_new_members} new {registry.__name__} records")
            save(
                new_members, ignore_conflicts=True
            )  # conflicts arising from existing records are ignored

        # create a new feature set from feature values using the same uid
        schema_self = Schema.from_values(member_uids, field=getattr(registry, field))
        if schema_self is None:
            if hasattr(registry, "organism_id"):
                logger.warning(
                    f"Schema is not transferred, check if organism is set correctly: {schema}"
                )
            continue
        # make sure the uid matches if schema is composed of same features
        if schema_self.hash == schema.hash:
            schema_self.uid = schema.uid
        logger.info(f"saving {slot} schema: {schema_self}")
        self._host.features._add_schema(schema_self, slot)


def make_external(self, feature: Feature) -> None:
    """Make a feature external, aka, remove feature from feature sets.

    Args:
        feature: `Feature` A feature record.

    """
    if not isinstance(feature, Feature):
        raise TypeError("feature must be a Feature record!")
    feature_sets = Schema.filter(features=feature).all()
    for fs in feature_sets:
        f = Feature.filter(uid=feature.uid).all()
        features_updated = fs.members.difference(f)
        if len(features_updated) > 0:
            # re-compute the hash of feature sets based on the updated members
            features_hash = hash_set({feature.uid for feature in features_updated})
            fs.hash = features_hash
            fs.n = len(features_updated)
            fs.save()
        # delete the link between the feature and the feature set
        Schema.features.through.objects.filter(
            feature_id=feature.id, schema_id=fs.id
        ).delete()
        # if no members are left in the schema, delete it
        if len(features_updated) == 0:
            logger.warning(f"deleting empty feature set: {fs}")
            fs.artifacts.set([])
            fs.delete()


@deprecated("_add_schema")
def add_schema(self, schema: Schema, slot: str) -> None:
    return self._add_schema(schema, slot)


@deprecated("_add_schema")
def add_feature_set(self, schema: Schema, slot: str) -> None:
    return self._add_schema(schema, slot)


@property
@deprecated("slots")
def _schema_by_slot(self):
    return self.slots


@property
def _feature_set_by_slot(self):
    return self.slots


# deprecated: feature set parsing


def parse_staged_feature_sets_from_anndata(
    adata: AnnData,
    var_field: FieldAttr | None = None,
    obs_field: FieldAttr = Feature.name,
    uns_field: FieldAttr | None = None,
    mute: bool = False,
    organism: str | SQLRecord | None = None,
) -> dict:
    data_parse = adata
    if not isinstance(adata, AnnData):  # is a path
        filepath = create_path(adata)  # returns Path for local
        if not isinstance(filepath, LocalPathClasses):
            from lamindb import settings
            from lamindb.core.storage._backed_access import backed_access

            using_key = settings._using_key
            data_parse = backed_access(filepath, using_key=using_key)
        else:
            data_parse = ad.read_h5ad(filepath, backed="r")
        type = "float"
    else:
        type = "float" if adata.X is None else serialize_pandas_dtype(adata.X.dtype)
    feature_sets = {}
    if var_field is not None:
        schema_var = Schema.from_values(
            data_parse.var.index,
            var_field,
            type=type,
            mute=mute,
            organism=organism,
            raise_validation_error=False,
        )
        if schema_var is not None:
            feature_sets["var"] = schema_var
    if obs_field is not None and len(data_parse.obs.columns) > 0:
        schema_obs = Schema.from_df(
            df=data_parse.obs,
            field=obs_field,
            mute=mute,
            organism=organism,
        )
        if schema_obs is not None:
            feature_sets["obs"] = schema_obs
    if uns_field is not None and len(data_parse.uns) > 0:
        validated_features = Feature.from_values(  # type: ignore
            data_parse.uns.keys(), field=uns_field, organism=organism
        )
        if len(validated_features) > 0:
            schema_uns = Schema(validated_features, dtype=None, otype="dict")
            feature_sets["uns"] = schema_uns
    return feature_sets


# no longer called from within curator
# might deprecate in the future?
def _add_set_from_df(
    self,
    field: FieldAttr = Feature.name,
    organism: str | None = None,
    mute: bool = False,
):
    """Add feature set corresponding to column names of DataFrame."""
    assert self._host.otype == "DataFrame"  # noqa: S101
    df = self._host.load(is_run_input=False)
    schema = Schema.from_df(
        df=df,
        field=field,
        mute=mute,
        organism=organism,
    )
    self._host._staged_feature_sets = {"columns": schema}
    self._host.save()


def _add_set_from_anndata(
    self,
    var_field: FieldAttr | None = None,
    obs_field: FieldAttr | None = Feature.name,
    uns_field: FieldAttr | None = None,
    mute: bool = False,
    organism: str | SQLRecord | None = None,
):
    """Add features from AnnData."""
    assert self._host.otype == "AnnData"  # noqa: S101

    # parse and register features
    adata = self._host.load(is_run_input=False)
    feature_sets = parse_staged_feature_sets_from_anndata(
        adata,
        var_field=var_field,
        obs_field=obs_field,
        uns_field=uns_field,
        mute=mute,
        organism=organism,
    )

    # link feature sets
    self._host._staged_feature_sets = feature_sets
    self._host.save()


def _add_set_from_mudata(
    self,
    var_fields: dict[str, FieldAttr] | None = None,
    obs_fields: dict[str, FieldAttr] | None = None,
    mute: bool = False,
    organism: str | SQLRecord | None = None,
):
    """Add features from MuData."""
    if obs_fields is None:
        obs_fields = {}
    assert self._host.otype == "MuData"  # noqa: S101

    # parse and register features
    mdata = self._host.load(is_run_input=False)
    feature_sets = {}

    obs_features = Feature.from_values(mdata.obs.columns)  # type: ignore
    if len(obs_features) > 0:
        feature_sets["obs"] = Schema(features=obs_features)
    for modality, field in var_fields.items():
        modality_fs = parse_staged_feature_sets_from_anndata(
            mdata[modality],
            var_field=field,
            obs_field=obs_fields.get(modality, Feature.name),
            mute=mute,
            organism=organism,
        )
        for k, v in modality_fs.items():
            feature_sets[f"['{modality}'].{k}"] = v

    # link feature sets
    self._host._staged_feature_sets = _unify_staged_feature_sets_by_hash(feature_sets)
    self._host.save()


def _add_set_from_spatialdata(
    self,
    sample_metadata_key: str,
    sample_metadata_field: FieldAttr = Feature.name,
    var_fields: dict[str, FieldAttr] | None = None,
    obs_fields: dict[str, FieldAttr] | None = None,
    mute: bool = False,
    organism: str | SQLRecord | None = None,
):
    """Add features from SpatialData."""
    obs_fields, var_fields = obs_fields or {}, var_fields or {}
    assert self._host.otype == "SpatialData"  # noqa: S101

    # parse and register features
    sdata = self._host.load(is_run_input=False)
    feature_sets = {}

    # sample features
    sample_features = Feature.from_values(
        sdata.get_attrs(key=sample_metadata_key, return_as="df", flatten=True).columns,
        field=sample_metadata_field,
    )  # type: ignore
    if len(sample_features) > 0:
        feature_sets[sample_metadata_key] = Schema(features=sample_features)

    # table features
    for table, field in var_fields.items():
        table_fs = parse_staged_feature_sets_from_anndata(
            sdata[table],
            var_field=field,
            obs_field=obs_fields.get(table, Feature.name),
            mute=mute,
            organism=organism,
        )
        for k, v in table_fs.items():
            feature_sets[f"['{table}'].{k}"] = v

    # link feature sets
    self._host._staged_feature_sets = _unify_staged_feature_sets_by_hash(feature_sets)
    self._host.save()


# mypy: ignore-errors
FeatureManager.__init__ = __init__
FeatureManager.__repr__ = __repr__
FeatureManager.describe = describe
FeatureManager.__getitem__ = __getitem__
FeatureManager.get_values = get_values
FeatureManager.slots = slots
FeatureManager.add_values = add_values_features
FeatureManager._add_schema = _add_schema
FeatureManager._accessor_by_registry = _accessor_by_registry
FeatureManager._add_from = _add_from
FeatureManager.filter = filter
FeatureManager.get = get
FeatureManager.make_external = make_external
FeatureManager.remove_values = remove_values

# deprecated
FeatureManager._add_set_from_df = _add_set_from_df
FeatureManager._add_set_from_anndata = _add_set_from_anndata
FeatureManager._add_set_from_mudata = _add_set_from_mudata
FeatureManager._add_set_from_spatialdata = _add_set_from_spatialdata
FeatureManager.add_schema = add_schema
FeatureManager.add_feature_set = add_feature_set
FeatureManager._schema_by_slot = _schema_by_slot
FeatureManager._feature_set_by_slot = _feature_set_by_slot
