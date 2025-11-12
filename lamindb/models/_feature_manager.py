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
from django.db.models import Aggregate, Subquery
from django.db.utils import IntegrityError
from lamin_utils import logger
from lamindb_setup.core.hashing import hash_set
from lamindb_setup.core.upath import create_path
from lamindb_setup.errors import ModuleWasntConfigured
from rich.table import Column, Table
from rich.text import Text
from rich.tree import Tree

from lamindb.core.storage import LocalPathClasses
from lamindb.errors import DoesNotExist, InvalidArgument, ValidationError
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
from ._django import get_artifact_or_run_with_related
from ._label_manager import _get_labels
from ._relations import (
    dict_related_model_to_related_name,
)
from .feature import Feature, FeatureValue, parse_dtype
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
    from lamindb.models.query_set import BasicQuerySet

    from .record import Record
    from .run import Run


def get_accessor_by_registry_(host: Artifact | Collection) -> dict:
    dictionary = {
        field.related_model.__get_name_with_module__(): field.name
        for field in host._meta.related_objects
    }
    dictionary["Feature"] = "features"
    dictionary["ULabel"] = "ulabels"
    dictionary["Record"] = "records"
    return dictionary


def get_schema_by_slot_(host: Artifact) -> dict[str, Schema]:
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
) -> BasicQuerySet:
    kwargs = {"artifact_id": host.id, "feature_id": feature.id}
    link_records = (
        getattr(host, host.features._accessor_by_registry[registry])  # type: ignore
        .through.objects.using(host._state.db)
        .filter(**kwargs)
    )
    return link_records


def get_schema_links(host: Artifact | Collection) -> BasicQuerySet:
    kwargs = {"artifact_id": host.id}
    links_schema = host.feature_sets.through.objects.filter(**kwargs)
    return links_schema


def get_link_attr(link: IsLink | type[IsLink], data: Artifact | Collection) -> str:
    link_model_name = link.__class__.__name__
    if link_model_name in {"Registry", "ModelBase"}:  # we passed the type of the link
        link_model_name = link.__name__  # type: ignore
    return link_model_name.replace(data.__class__.__name__, "").lower()


def strip_cat(feature_dtype: str) -> str:
    if "cat[" in feature_dtype:
        parts = feature_dtype.split("cat[")
        dtype_stripped_cat = "".join(
            part[:-1] if i != 0 else part for i, part in enumerate(parts)
        )
    else:
        dtype_stripped_cat = feature_dtype
    return dtype_stripped_cat


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
    if related_data is None:
        if self.__class__.__name__ in {"Artifact", "Run", "Record"}:
            artifact_meta = get_artifact_or_run_with_related(
                self, include_feature_link=True, include_m2m=True
            )
            related_data = artifact_meta.get("related_data", {})
        else:
            related_data = {}

    # Process m2m data
    m2m_data = related_data.get("m2m", {}) if related_data else {}
    # e.g. m2m_data = {'tissues': {1: {'id': 1, 'uid': '1fIFAQJY', 'abbr': None, 'name': 'brain', 'tissue': 1, 'feature': 1, 'ontology_id': 'UBERON:0000955', 'tissue_display': 'brain'}, 10: {'id': 2, 'uid': '7Tt4iEKc', 'abbr': None, 'name': 'lung', 'tissue': 10, 'feature': 1, 'ontology_id': 'UBERON:0002048', 'tissue_display': 'lung'}}, 'cell_types': {1: {'id': 1, 'uid': '3QnZfoBk', 'abbr': None, 'name': 'neuron', 'feature': 2, 'celltype': 1, 'ontology_id': 'CL:0000540', 'celltype_display': 'neuron'}}}
    # e.g. {'tissue': {1: {'id': 1, 'uid': '1fIFAQJY', 'abbr': None, 'name': 'brain', 'tissue': 1, 'feature': 1, 'ontology_id': 'UBERON:0000955', 'tissue_display': 'brain'}, 10: {'id': 2, 'uid': '7Tt4iEKc', 'abbr': None, 'name': 'lung', 'tissue': 10, 'feature': 1, 'ontology_id': 'UBERON:0002048', 'tissue_display': 'lung'}}, 'celltype': {1: {'id': 1, 'uid': '3QnZfoBk', 'abbr': None, 'name': 'neuron', 'feature': 2, 'celltype': 1, 'ontology_id': 'CL:0000540', 'celltype_display': 'neuron'}}}
    # integers are the ids of the related labels
    m2m_name = {}
    if not self.__class__.__name__ == "Record":
        for related_name, values in m2m_data.items():
            link_model = getattr(self.__class__, related_name).through
            related_model_name = link_model.__name__.replace(
                self.__class__.__name__, ""
            ).lower()
            m2m_name[related_model_name] = values
    else:
        m2m_name = related_data.get("m2m", {})

    # Get feature information
    links_data = related_data.get("link", {}) if related_data else {}
    # e.g. feature_dict = {1: ('tissue', 'cat[bionty.Tissue.ontology_id]'), 2: ('cell_type', 'cat[bionty.CellType]')}
    feature_dict = {
        id: (name, dtype)
        for id, name, dtype in Feature.objects.using(self._state.db).values_list(
            "id", "name", "dtype"
        )
    }

    # Build result dictionary
    result = {}  # type: ignore
    for link_name, link_values in links_data.items():
        related_name = link_name.removeprefix("links_").replace("_", "")
        if not link_values:
            continue
        # sort by the order on the link table, important for list dtypes
        for link_value in sorted(link_values, key=lambda x: x.get("id")):
            feature_id = link_value.get("feature")
            if feature_id is None:
                continue
            feature_name, feature_dtype = feature_dict.get(feature_id)
            feature_field = parse_dtype(feature_dtype)[0]["field_str"]
            if not self.__class__.__name__ == "Record":
                label_id = link_value.get(related_name)
                label_name = (
                    m2m_name.get(related_name, {}).get(label_id).get(feature_field)
                )
            else:
                label_name = link_value.get(feature_field)
            if label_name:
                dict_key = (feature_name, feature_dtype)
                if dict_key not in result:
                    result[dict_key] = (
                        set() if not feature_dtype.startswith("list[cat") else []
                    )
                if feature_dtype.startswith("list[cat"):
                    result[dict_key].append(label_name)
                else:
                    result[dict_key].add(label_name)
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
                feature_field = parse_dtype(feature.dtype)[0]["field_str"]
                link_attr = get_link_attr(link, self)
                label = getattr(link, link_attr)
                label_name = getattr(label, feature_field)
                result[(feature.name, feature.dtype)].add(label_name)

    return dict(result)


def _get_non_categoricals(
    self,
) -> dict[tuple[str, str], set[Any]]:
    """Get non-categorical features and their values."""
    from .artifact import Artifact
    from .record import Record
    from .run import Run

    non_categoricals = {}

    if self.id is not None and isinstance(self, (Artifact, Run, Record)):
        if isinstance(self, Record):
            _feature_values = self.values_json.values(
                "feature__name", "feature__dtype", "value"
            ).order_by("feature__name")
        else:
            _feature_values = (
                self._feature_values.values("feature__name", "feature__dtype")
                .annotate(values=custom_aggregate("value", self._state.db))
                .order_by("feature__name")
            )

        for fv in _feature_values:
            feature_name = fv["feature__name"]
            feature_dtype = fv["feature__dtype"]
            if isinstance(self, Record):
                values = fv["value"]
            else:
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


def get_features_data(
    self: Artifact | Run | Record,
    related_data: dict | None = None,
    to_dict: bool = False,
    external_only: bool = False,
):
    from .artifact import Artifact

    dictionary: dict[str, Any] = {}

    if self._state.adding:
        if to_dict:
            return dictionary
        else:
            raise NotImplementedError

    # feature sets
    schema_data: dict[str, tuple[str, list[str]]] = {}
    feature_data: dict[str, tuple[str, list[str]]] = {}
    if not to_dict and isinstance(self, Artifact):
        if self.id is not None and connections[self._state.db].vendor == "postgresql":
            if not related_data:
                artifact_meta = get_artifact_or_run_with_related(
                    self,
                    include_schema=True,
                    include_m2m=True,
                    include_feature_link=True,
                )
                related_data = artifact_meta.get("related_data", {})
            fs_data = related_data.get("schemas", {}) if related_data else {}
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
                    slot: (schema, schema.n)  # type: ignore
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

    internal_feature_names = {}
    if isinstance(self, Artifact):
        inferred_schemas = self.feature_sets.filter(itype="Feature")
        if len(inferred_schemas) > 0:
            for schema in inferred_schemas:
                internal_feature_names.update(
                    dict(schema.members.values_list("name", "dtype"))
                )

    # categorical feature values
    # Get the categorical data using the appropriate method
    # e.g. categoricals = {('tissue', 'cat[bionty.Tissue.ontology_id]'): {'brain'}, ('cell_type', 'cat[bionty.CellType]'): {'neuron'}}
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

    internal_feature_labels = {}
    external_data = []
    for features, is_list_type in [(categoricals, False), (non_categoricals, True)]:
        for (feature_name, feature_dtype), values in sorted(features.items()):
            # Handle dictionary conversion
            if to_dict:
                if feature_dtype.startswith("list[cat"):
                    dict_value = values  # is already a list
                else:
                    dict_value = values if len(values) > 1 else next(iter(values))
                dictionary[feature_name] = dict_value
                continue

            # Format message
            printed_values = (
                _format_values(sorted(values), n=10, quotes=False)
                if not is_list_type or not feature_dtype.startswith("list")
                else str(values)  # need to convert to string
            )

            # Sort into internal/external
            feature_info = (
                feature_name,
                Text(strip_cat(feature_dtype), style="dim"),
                printed_values,
            )
            if feature_name in internal_feature_names:
                internal_feature_labels[feature_name] = feature_info
            else:
                external_data.append(feature_info)

    if to_dict:
        if external_only:
            return {
                k: v for k, v in dictionary.items() if k not in internal_feature_names
            }
        else:
            return dictionary
    else:
        return (
            internal_feature_labels,
            feature_data,
            schema_data,
            internal_feature_names,
            external_data,
        )


def describe_features(
    self: Artifact | Run | Record,
    related_data: dict | None = None,
) -> tuple[Tree | None, Tree | None]:
    """Describe features of an artifact or collection."""
    if self._state.adding:
        return None, None
    (
        internal_feature_labels,
        feature_data,
        schema_data,
        internal_feature_names,
        external_data,
    ) = get_features_data(
        self,
        related_data=related_data,
    )

    # Dataset features section
    # internal features that contain labels (only `Feature` features contain labels)
    internal_feature_labels_slot: dict[str, list] = {}
    for feature_name, feature_row in internal_feature_labels.items():
        slot, _ = feature_data.get(feature_name)
        internal_feature_labels_slot.setdefault(slot, []).append(feature_row)

    dataset_features_tree_children = []
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
                            strip_cat(internal_feature_names.get(feature_name)),
                            style="dim",
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
                            strip_cat(
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
        schema_itype = f" {schema.itype}" if schema.itype != "Feature" else ""
        dataset_features_tree_children.append(
            _create_feature_table(
                Text.assemble(
                    (slot, "violet"),
                    (f" ({schema.n}{schema_itype})", "dim"),
                ),
                "",
                feature_rows,
                show_header=True,
            )
        )
    # external features
    external_features_tree_children = []
    if external_data:
        external_features_tree_children.append(
            _create_feature_table(
                "",
                "",
                external_data,
            )
        )

    # trees
    dataset_features_tree = None
    if dataset_features_tree_children:
        dataset_features_tree = Tree(
            Text("Dataset features", style="bold bright_magenta")
        )
        for child in dataset_features_tree_children:
            dataset_features_tree.add(child)
    external_features_tree = None
    if external_features_tree_children:
        external_features_text = (
            "External features"
            if (
                self.__class__.__name__ == "Artifact" and dataset_features_tree_children
            )
            else "Features"
        )
        external_features_tree = Tree(
            Text(external_features_text, style="bold dark_orange")
        )
        for child in external_features_tree_children:
            external_features_tree.add(child)
    return dataset_features_tree, external_features_tree


def infer_feature_type_convert_json(
    key: str, value: Any, mute: bool = False
) -> tuple[str, Any, str]:
    from lamindb.base.dtypes import is_valid_datetime_str

    message = ""
    if isinstance(value, bool):
        return "bool", value, message
    elif isinstance(value, int):
        return "int", value, message
    elif isinstance(value, float):
        return "float", value, message
    elif isinstance(value, datetime):
        return "datetime", value.isoformat(), message
    elif isinstance(value, date):
        return "date", value.isoformat(), message
    elif isinstance(value, str):
        if datetime_str := is_valid_datetime_str(value):
            dt_type = (
                "date" if len(value) == 10 else "datetime"
            )  # YYYY-MM-DD is exactly 10 characters
            sanitized_value = datetime_str[:10] if dt_type == "date" else datetime_str  # type: ignore
            return dt_type, sanitized_value, message  # type: ignore
        else:
            return "cat ? str", value, message
    elif isinstance(value, SQLRecord):
        return (f"cat[{value.__class__.__get_name_with_module__()}]", value, message)
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
            first_element = next(iter(value))
            first_element_type = type(first_element)
            # check that all elements are of the same type
            if all(isinstance(elem, first_element_type) for elem in value):
                if first_element_type is bool:
                    return "list[bool]", value, message
                elif first_element_type is int:
                    return "list[int]", value, message
                elif first_element_type is float:
                    return "list[float]", value, message
                elif first_element_type is str:
                    return ("list[cat ? str]", value, message)
                elif isinstance(first_element, SQLRecord):
                    return (
                        f"list[cat[{first_element_type.__get_name_with_module__()}]]",
                        value,
                        message,
                    )
    if not mute:
        logger.warning(f"cannot infer feature type of: {value}, returning '?")
    return "?", value, message


def filter_base(
    queryset: BasicQuerySet,
    _skip_validation: bool = True,
    **expression,
) -> BasicQuerySet:
    from lamindb.models import Artifact, BasicQuerySet, QuerySet

    # not QuerySet but only BasicQuerySet
    assert isinstance(queryset, BasicQuerySet) and not isinstance(queryset, QuerySet)  # noqa: S101

    registry = queryset.model
    db = queryset.db

    model = Feature
    value_model = FeatureValue
    keys_normalized = [key.split("__")[0] for key in expression]
    if not _skip_validation:
        validated = model.connect(db).validate(keys_normalized, field="name", mute=True)
        if sum(validated) != len(keys_normalized):
            raise ValidationError(
                f"Some keys in the filter expression are not registered as features: {np.array(keys_normalized)[~validated]}"
            )
    new_expression = {}
    features = model.connect(db).filter(name__in=keys_normalized).distinct()
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
                if registry is Artifact:
                    from .artifact import ArtifactFeatureValue

                    if value:  # True
                        return queryset.exclude(
                            id__in=Subquery(
                                ArtifactFeatureValue.objects.filter(
                                    featurevalue__feature=feature
                                ).values("artifact_id")
                            )
                        )
                    else:
                        return queryset.exclude(
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
                if registry is Artifact:
                    result = parse_dtype(feature.dtype)[0]
                    kwargs = {
                        f"links_{result['registry'].__name__.lower()}__feature": feature
                    }
                    if value:  # True
                        return queryset.exclude(**kwargs)
                    else:
                        return queryset.filter(**kwargs)
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
                    labels = result["registry"].connect(db).filter(**expression)
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
    if not new_expression:
        raise NotImplementedError
    return queryset.filter(**new_expression)


def filter_with_features(
    queryset: BasicQuerySet, *queries, **expressions
) -> BasicQuerySet:
    from lamindb.models import Artifact, BasicQuerySet, QuerySet

    if isinstance(queryset, QuerySet):
        # need to avoid infinite recursion because
        # filter_with_features is called in queryset.filter otherwise
        filter_kwargs = {"_skip_filter_with_features": True}
    else:
        filter_kwargs = {}

    registry = queryset.model

    if registry is Artifact and not any(e.startswith("kind") for e in expressions):
        exclude_kwargs = {"kind": "__lamindb_run__"}
    else:
        exclude_kwargs = {}

    if expressions:
        keys_normalized = [key.split("__")[0] for key in expressions]
        field_or_feature = keys_normalized[0]
        if field_or_feature in registry.__get_available_fields__():
            qs = queryset.filter(*queries, **expressions, **filter_kwargs)
        elif all(
            features_validated := Feature.objects.using(queryset.db).validate(
                keys_normalized, field="name", mute=True
            )
        ):
            # filter_base requires qs to be BasicQuerySet
            qs = filter_base(
                queryset._to_class(BasicQuerySet, copy=True),
                _skip_validation=True,
                **expressions,
            )._to_class(type(queryset), copy=False)
            qs = qs.filter(*queries, **filter_kwargs)
        else:
            features = ", ".join(sorted(np.array(keys_normalized)[~features_validated]))
            message = f"feature names: {features}"
            avail_fields = registry.__get_available_fields__()
            if "_branch_code" in avail_fields:
                avail_fields.remove("_branch_code")  # backward compat
            fields = ", ".join(sorted(avail_fields))
            raise InvalidArgument(
                f"You can query either by available fields: {fields}\n"
                f"Or fix invalid {message}"
            )
    else:
        qs = queryset.filter(*queries, **filter_kwargs)

    return qs.exclude(**exclude_kwargs) if exclude_kwargs else qs


# for deprecated functionality
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


# for deprecated functionality
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
        dtype = "float"
    else:
        dtype = "float" if adata.X is None else serialize_pandas_dtype(adata.X.dtype)
    feature_sets = {}
    if var_field is not None:
        schema_var = Schema.from_values(
            data_parse.var.index,
            var_field,
            dtype=dtype,
            mute=mute,
            organism=organism,
            raise_validation_error=False,
        )
        if schema_var is not None:
            feature_sets["var"] = schema_var
    if obs_field is not None and len(data_parse.obs.columns) > 0:
        schema_obs = Schema.from_dataframe(
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


class FeatureManager:
    """Feature manager."""

    def __init__(self, host: Artifact | Run | Record):
        self._host = host
        self._slots: dict[str, Schema] | None = None
        self._accessor_by_registry_ = None

    def __repr__(self) -> str:
        return self.describe(return_str=True)  # type: ignore

    def describe(self, return_str: bool = False) -> str | None:
        """Pretty print features.

        This is what `artifact.describe()` calls under the hood.
        """
        dataset_features_tree, external_features_tree = describe_features(self._host)  # type: ignore
        tree = describe_header(self._host)
        if dataset_features_tree:
            tree.add(dataset_features_tree)
        if external_features_tree:
            tree.add(external_features_tree)
        return format_rich_tree(
            tree, fallback="no linked features", return_str=return_str
        )

    def get_values(self) -> dict[str, Any]:
        """Get features as a dictionary."""
        return get_features_data(self._host, to_dict=True)  # type: ignore

    @deprecated("slots[slot].members")
    def __getitem__(self, slot) -> BasicQuerySet:
        if slot not in self.slots:
            raise ValueError(
                f"No linked feature set for slot: {slot}\nDid you get validation"
                " warnings? Only features that match registered features get validated"
                " and linked."
            )
        schema = self.slots[slot]
        orm_name = schema.itype
        return getattr(schema, self._accessor_by_registry[orm_name]).all()

    @property
    def slots(self) -> dict[str, Schema]:
        """Features by schema slot.

        Example::

            artifact.features.slots
            #> {'var': <Schema: var>, 'obs': <Schema: obs>}
        """
        if self._slots is None:
            self._slots = get_schema_by_slot_(self._host)
        return self._slots

    @property
    def _accessor_by_registry(self):
        """Accessor by registry."""
        if self._accessor_by_registry_ is None:
            self._accessor_by_registry_ = get_accessor_by_registry_(self._host)
        return self._accessor_by_registry_

    def _add_label_feature_links(
        self,
        features_labels,
        *,
        label_ref_is_name: bool | None = None,
        feature_ref_is_name: bool | None = None,
    ):
        host_name = self._host.__class__.__name__.lower()
        host_is_record = host_name == "record"
        related_names = dict_related_model_to_related_name(self._host.__class__)
        if host_is_record:
            # the related model to related name is ambiguous if two M2M exist
            related_names["Record"] = "components"
            related_names["Project"] = "linked_projects"
            related_names["Artifact"] = "linked_artifacts"
            related_names["Run"] = "linked_runs"
        for class_name, registry_features_labels in features_labels.items():
            related_name = related_names[class_name]  # e.g., "ulabels"
            IsLink = getattr(self._host, related_name).through
            if host_is_record:
                field_name = "value_id"
            else:
                field_name = (
                    f"{get_link_attr(IsLink, self._host)}_id"  # e.g., ulabel_id
                )
            if host_name == "artifact":
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
            else:  # Run
                links = [
                    IsLink(
                        **{
                            f"{host_name}_id": self._host.id,
                            "feature_id": feature.id,
                            field_name: label.id,
                        }
                    )
                    for (feature, label) in registry_features_labels
                ]
            # a link might already exist
            try:
                save(links, ignore_conflicts=False)
            except Exception:
                save(links, ignore_conflicts=True)

    def _get_feature_records(self, dictionary, feature_field):
        from .._tracked import get_current_tracked_run

        registry = feature_field.field.model
        keys = list(dictionary.keys())
        feature_records = registry.from_values(keys, field=feature_field, mute=True)
        if len(feature_records) != len(keys):
            not_validated_keys = [
                key for key in keys if key not in feature_records.to_list("name")
            ]
            not_validated_keys_dtype_message = [
                (key, infer_feature_type_convert_json(key, dictionary[key]))
                for key in not_validated_keys
            ]
            run = get_current_tracked_run()
            if run is not None:
                name = f"{run.transform.type}[{run.transform.key}]"
                type_hint = f"""  feature_type = ln.Feature(name='{name}', is_type=True).save()"""
                elements = [type_hint]
                type_kwarg = ", type=feature_type"
            else:
                elements = []
                type_kwarg = ""
            elements += [
                f"  ln.Feature(name='{key}', dtype='{dtype}'{type_kwarg}).save(){message}"
                for key, (dtype, _, message) in not_validated_keys_dtype_message
            ]
            hint = "\n".join(elements)
            msg = (
                f"These keys could not be validated: {not_validated_keys}\n"
                f"Here is how to create a feature:\n\n{hint}"
            )
            raise ValidationError(msg)
        return feature_records

    def add_values(
        self,
        values: dict[str, str | int | float | bool],
        feature_field: FieldAttr = Feature.name,
        schema: Schema = None,
    ) -> None:
        """Add values for features.

        Args:
            values: A dictionary of keys (features) & values (labels, strings, numbers, booleans, datetimes, etc.).
                If a value is `None`, it will be skipped.
            feature_field: The field of a registry to map the keys of the `values` dictionary.
            schema: Schema to validate against.
        """
        from lamindb.curators.core import ExperimentalDictCurator

        host_is_record = self._host.__class__.__name__ == "Record"
        host_is_artifact = self._host.__class__.__name__ == "Artifact"
        # rename to distinguish from the values inside the dict
        dictionary = values
        keys = dictionary.keys()
        if isinstance(keys, DICT_KEYS_TYPE):
            keys = list(keys)  # type: ignore
        # deal with other cases later
        assert all(isinstance(key, str) for key in keys)  # noqa: S101
        if (
            host_is_record
            and self._host.type is not None
            and self._host.type.schema is not None  # type: ignore
        ):
            assert schema is None, "Cannot pass schema if record.type has schema."
            schema = self._host.type.schema  # type: ignore
        if host_is_artifact:
            if self._get_external_schema():
                raise ValueError("Cannot add values if artifact has external schema.")
        if schema is not None:
            ExperimentalDictCurator(values, schema).validate()
            feature_records = schema.members.filter(name__in=keys)
        else:
            feature_records = self._get_feature_records(dictionary, feature_field)
        return self._add_values(feature_records, dictionary)

    def _add_values(self, feature_records, dictionary):
        from ..base.dtypes import is_iterable_of_sqlrecord
        from .can_curate import CanCurate
        from .record import RecordJson

        host_is_record = self._host.__class__.__name__ == "Record"
        features_labels = defaultdict(list)
        feature_json_values = []
        not_validated_values: dict[str, tuple[str, list[str]]] = {}
        for feature in feature_records:
            value = dictionary[feature.name]
            if value is None:
                continue
            inferred_type, converted_value, _ = infer_feature_type_convert_json(
                feature.name,
                value,
                mute=True,
            )
            if feature.dtype == "num":
                if inferred_type not in {"int", "float"}:
                    raise TypeError(
                        f"Value for feature '{feature.name}' with dtype {feature.dtype} must be a number, but is {value} with dtype {inferred_type}"
                    )
            elif feature.dtype.startswith("cat"):
                if inferred_type != "?":
                    if not (
                        inferred_type.startswith("cat")
                        or inferred_type == "list[cat ? str]"
                        or isinstance(value, SQLRecord)
                        or is_iterable_of_sqlrecord(value)
                    ):
                        raise TypeError(
                            f"Value for feature '{feature.name}' with dtype '{feature.dtype}' must be a string or record, but is {value} with dtype {inferred_type}"
                        )
            elif (
                (feature.dtype == "str" and inferred_type != "cat ? str")
                or (feature.dtype == "list[str]" and inferred_type != "list[cat ? str]")
                or (
                    feature.dtype.startswith("list[cat")
                    and not inferred_type.startswith("list[cat")
                )
                or (
                    feature.dtype not in {"str", "list[str]"}
                    and not feature.dtype.startswith("list[cat")
                    and feature.dtype != inferred_type
                )
            ):
                raise ValidationError(
                    f"Expected dtype for '{feature.name}' is '{feature.dtype}', got '{inferred_type}'"
                )
            if not (
                feature.dtype.startswith("cat") or feature.dtype.startswith("list[cat")
            ):
                filter_kwargs = {"feature": feature, "value": converted_value}
                if host_is_record:
                    filter_kwargs["record"] = self._host
                    feature_value = RecordJson(**filter_kwargs)
                else:
                    feature_value, _ = FeatureValue.get_or_create(**filter_kwargs)
                feature_json_values.append(feature_value)
            else:
                if isinstance(value, SQLRecord) or is_iterable_of_sqlrecord(value):
                    if isinstance(value, SQLRecord):
                        label_records = [value]
                    else:
                        label_records = value  # type: ignore
                    for record in label_records:
                        if record._state.adding:
                            raise ValidationError(
                                f"Please save {record} before annotation."
                            )
                        features_labels[
                            record.__class__.__get_name_with_module__()
                        ].append((feature, record))
                else:
                    if isinstance(value, str):
                        values = [value]  # type: ignore
                    else:
                        values = value  # type: ignore
                    if feature.dtype == "cat":
                        feature.dtype += "[ULabel]"
                        feature.save()
                        result = {
                            "registry_str": "ULabel",
                            "registry": ULabel,
                            "field": ULabel.name,
                        }
                    else:
                        result = parse_dtype(feature.dtype)[0]
                    if issubclass(result["registry"], CanCurate):  # type: ignore
                        validated = result["registry"].validate(  # type: ignore
                            values, field=result["field"], mute=True
                        )
                        values_array = np.array(values)
                        validated_values = values_array[validated]
                        key = result["registry_str"]
                        if validated.sum() != len(values):
                            not_validated_values[result["registry_str"]] = (  # type: ignore
                                result["field_str"],
                                values_array[~validated].tolist(),
                            )
                        label_records = result["registry"].from_values(  # type: ignore
                            validated_values, field=result["field"], mute=True
                        )
                    else:
                        label_records = result["registry"].filter(  # type: ignore
                            **{f"{result['field_str']}__in": values}
                        )
                        if len(label_records) != len(values):
                            raise ValidationError(
                                f"Some of these values for {result['registry_str']} do not exist: {values}"
                            )
                    features_labels[result["registry_str"]] += [  # type: ignore
                        (feature, label_record) for label_record in label_records
                    ]
        if not_validated_values:
            hint = ""
            for key, (field, values_list) in not_validated_values.items():
                key_str = "ln.Record" if key == "Record" else key
                create_true = ", create=True" if "bionty." not in key else ""
                hint += f"  records = {key_str}.from_values({values_list}, field='{field}'{create_true}).save()\n"
            msg = (
                f"These values could not be validated: {dict(not_validated_values)}\n"
                f"Here is how to create records for them:\n\n{hint}"
            )
            raise ValidationError(msg)
        if features_labels:
            self._add_label_feature_links(features_labels)
        if feature_json_values and host_is_record:
            save(feature_json_values)
        elif feature_json_values:
            to_insert_feature_values = [
                record for record in feature_json_values if record._state.adding
            ]
            if to_insert_feature_values:
                save(to_insert_feature_values)
            links = [
                self._host._feature_values.through(
                    **{
                        f"{self._host.__class__.__name__.lower()}_id": self._host.id,
                        "featurevalue_id": json_value.id,
                    }
                )
                for json_value in feature_json_values
            ]
            # a link might already exist, hence ignore_conflicts is needed
            save(links, ignore_conflicts=True)

    def set_values(
        self,
        values: dict[str, str | int | float | bool],
        feature_field: FieldAttr = Feature.name,
        schema: Schema = None,
    ) -> None:
        """Set values for features.

        Like `add_values`, but first removes all existing external feature annotations.

        Args:
            values: A dictionary of keys (features) & values (labels, strings, numbers, booleans, datetimes, etc.).
                If a value is `None`, it will be skipped.
            feature_field: The field of a registry to map the keys of the `values` dictionary.
            schema: Schema to validate against.
        """
        from lamindb.curators.core import ExperimentalDictCurator

        host_is_record = self._host.__class__.__name__ == "Record"
        host_is_artifact = self._host.__class__.__name__ == "Artifact"
        # rename to distinguish from the values inside the dict
        dictionary = values
        keys = dictionary.keys()
        if isinstance(keys, DICT_KEYS_TYPE):
            keys = list(keys)  # type: ignore
        # deal with other cases later
        assert all(isinstance(key, str) for key in keys)  # noqa: S101
        if (
            host_is_record
            and self._host.type is not None
            and self._host.type.schema is not None  # type: ignore
        ):
            assert schema is None, "Cannot pass schema if record.type has schema."
            schema = self._host.type.schema  # type: ignore
        if host_is_artifact:
            schema = self._get_external_schema()
        if schema is not None:
            ExperimentalDictCurator(values, schema).validate()
            feature_records = schema.members.filter(name__in=keys)
        else:
            feature_records = self._get_feature_records(dictionary, feature_field)
        self._remove_values()
        self._add_values(
            feature_records,
            dictionary=dictionary,
        )

    def _get_external_schema(self) -> Schema | None:
        external_schema = None
        if self._host.otype is None:
            external_schema = self._host.schema
        elif self._host.schema is not None:
            external_schema = self._host.schema.slots.get("__external__", None)
        return external_schema

    def remove_values(
        self,
        feature: str | Feature | list[str | Feature] = None,
        *,
        value: Any | None = None,
    ) -> None:
        """Remove values for features.

        Args:
            feature: Indicate one or several features for which to remove values.
                If `None`, values for all external features will be removed.
            value: An optional value to restrict removal to a single value.
        """
        host_name = self._host.__class__.__name__.lower()
        host_is_artifact = host_name == "artifact"

        if host_is_artifact:
            external_schema = self._get_external_schema()
            if external_schema is not None:
                raise ValueError(
                    "Cannot remove values if artifact has external schema."
                )
        return self._remove_values(
            feature,
            value=value,
        )

    def _remove_values(
        self,
        feature: str | Feature | list[str | Feature] = None,
        *,
        value: Any | None = None,
    ) -> None:
        from django.apps import apps

        host_name = self._host.__class__.__name__.lower()
        host_is_record = host_name == "record"
        host_is_artifact = host_name == "artifact"

        if feature is None:
            features = get_features_data(
                self._host, to_dict=True, external_only=True
            ).keys()
        elif not isinstance(feature, list):
            features = [feature]
        else:
            features = feature
        for feature in features:
            if isinstance(feature, str):
                feature_record = Feature.get(name=feature)
            else:
                feature_record = feature
            if host_is_artifact:
                for schema in self.slots.values():
                    if feature_record in schema.members:
                        raise ValueError("Cannot remove values for dataset features.")
            filter_kwargs = {"feature": feature_record}
            none_message = f"with value {value!r} " if value is not None else ""
            if feature_record.dtype.startswith(("cat[", "list[cat")):  # type: ignore
                feature_registry = parse_dtype(feature_record.dtype)[0]["registry_str"]
                if "." in feature_registry:
                    parts = feature_registry.split(".")
                    app_label = parts[0]
                    entity_name = parts[-1]
                else:
                    app_label = "lamindb"
                    entity_name = feature_registry
                host_name = self._host.__class__.__name__
                link_model_name = f"{host_name}{entity_name}"
                link_model = apps.get_model(app_label, link_model_name)
                filter_kwargs[host_name.lower()] = self._host
                if value is not None:
                    if not isinstance(value, SQLRecord):
                        raise TypeError(
                            f"Expected a record for removing categorical feature value, "
                            f"got {value} of type {type(value)}"
                        )
                    assert not host_is_record, "Only artifacts support passing a value."
                    filter_kwargs[entity_name.lower()] = value
                link_records = link_model.objects.filter(**filter_kwargs)
                if not link_records.exists():
                    value_msg = f"with value {value!r} " if value is not None else ""
                    logger.warning(
                        f"no feature '{feature_record.name}' {value_msg}found on "
                        f"{host_name.lower()} '{self._host.uid}'!"
                    )
                    return
                link_records.delete()
            else:
                if value is not None:
                    filter_kwargs["value"] = value
                if host_is_record:
                    feature_values = self._host.values_json.filter(**filter_kwargs)
                else:
                    feature_values = self._host._feature_values.filter(**filter_kwargs)
                if not feature_values.exists():
                    logger.warning(
                        f"no feature '{feature_record.name}' {none_message}found on {self._host.__class__.__name__.lower()} '{self._host.uid}'!"
                    )
                    return
                if host_is_record:
                    feature_values.delete(permanent=True)
                else:
                    # the below might leave a dangling feature_value record
                    # but we don't want to pay the price of making another query just to remove this annotation
                    # we can clean the FeatureValue registry periodically if we want to
                    self._host._feature_values.remove(*feature_values)

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

    def _add_from(self, data: Artifact | Collection, transfer_logs: dict = None):
        """Transfer features from a artifact or collection."""
        # This only covers feature sets
        if transfer_logs is None:
            transfer_logs = {"mapped": [], "transferred": [], "run": None}
        from lamindb import settings

        using_key = settings._using_key
        for slot, schema in data.features.slots.items():  # type: ignore
            try:
                members = schema.members
            except ModuleWasntConfigured as err:
                logger.warning(f"skipping transfer of {slot} schema because {err}")
                continue
            if len(members) == 0:
                continue
            if len(members) > settings.annotation.n_max_records:
                logger.warning(
                    f"skipping creating {len(members)} > {settings.annotation.n_max_records} new {members[0].__class__.__name__} records"
                )
                schema_self = schema
                schema_exists = Schema.filter(hash=schema_self.hash).one_or_none()
                if schema_exists is not None:
                    schema_self = schema_exists
                else:
                    schema_self.save()
            else:
                registry = members[0].__class__
                # note here the features are transferred based on an unique field
                field = REGISTRY_UNIQUE_FIELD.get(registry.__name__.lower(), "uid")
                # this will be e.g. be a list of ontology_ids or uids
                member_uids = list(members.values_list(field, flat=True))
                validated = registry.validate(member_uids, field=field, mute=True)
                new_members_uids = list(compress(member_uids, ~validated))
                new_members = members.filter(**{f"{field}__in": new_members_uids})
                n_new_members = len(new_members)
                if len(members) > settings.annotation.n_max_records:
                    logger.warning(
                        f"skipping creating {n_new_members} > {settings.annotation.n_max_records} new {registry.__name__} records"
                    )
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
                            feature,
                            using_key,
                            transfer_fk=False,
                            transfer_logs=transfer_logs,
                        )
                    save(
                        new_members, ignore_conflicts=True
                    )  # conflicts arising from existing records are ignored

                # create a new feature set from feature values using the same uid
                schema_self = Schema.from_values(
                    member_uids, field=getattr(registry, field)
                )
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
            try:
                self._host.features._add_schema(schema_self, slot)
            except IntegrityError:
                logger.warning(
                    f"updating annotation of artifact {self._host.uid} with feature set for slot: {slot}"
                )
                self._host.feature_sets.through.objects.get(
                    artifact_id=self._host.id, slot=slot
                ).delete()
                self._host.features._add_schema(schema_self, slot)

    def make_external(self, feature: Feature) -> None:
        """Make a feature external.

        This removes a feature from `artifact.feature_sets` and thereby no longer marks it
        as a dataset feature.

        Args:
            feature: A feature.
        """
        if not isinstance(feature, Feature):
            raise TypeError("feature must be a Feature record!")
        feature_sets = Schema.filter(features=feature)
        for fs in feature_sets:
            f = Feature.filter(uid=feature.uid)
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

    # no longer called from within curator
    # deprecated
    def _add_set_from_df(
        self,
        field: FieldAttr = Feature.name,
        organism: str | None = None,
        mute: bool = False,
    ):
        """Add feature set corresponding to column names of DataFrame."""
        assert self._host.otype == "DataFrame"  # noqa: S101
        df = self._host.load(is_run_input=False)
        schema = Schema.from_dataframe(
            df=df,
            field=field,
            mute=mute,
            organism=organism,
        )
        self._host._staged_feature_sets = {"columns": schema}
        self._host.save()

    # deprecated
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

    # deprecated
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
        self._host._staged_feature_sets = _unify_staged_feature_sets_by_hash(
            feature_sets
        )
        self._host.save()

    # deprecated
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
            sdata.get_attrs(
                key=sample_metadata_key, return_as="df", flatten=True
            ).columns,
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
        self._host._staged_feature_sets = _unify_staged_feature_sets_by_hash(
            feature_sets
        )
        self._host.save()
