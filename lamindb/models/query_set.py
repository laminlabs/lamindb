"""Models library.

...

Query sets & managers
---------------------
.. autoclass:: BasicQuerySet
.. autoclass:: QuerySet
.. autoclass:: ArtifactSet
.. autoclass:: QueryManager
.. autoclass:: lamindb.models.query_set.BiontyDB
.. autoclass:: lamindb.models.query_set.PertdbDB

...
"""

from __future__ import annotations

import ast
import re
import warnings
from collections import UserList, defaultdict
from collections.abc import Iterable
from collections.abc import Iterable as IterableType
from importlib import import_module
from typing import TYPE_CHECKING, Any, Generic, NamedTuple, TypeVar, final

import lamindb_setup as ln_setup
import pandas as pd
from django.core.exceptions import FieldError
from django.db import models
from django.db.models import (
    F,
    FilteredRelation,
    ForeignKey,
    ManyToManyField,
    Q,
    Subquery,
)
from django.db.models.fields.related import ForeignObjectRel
from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core import deprecated
from lamindb_setup.core._docs import doc_args

from ..errors import DoesNotExist, MultipleResultsFound
from ._is_versioned import IsVersioned
from .can_curate import CanCurate, _inspect, _standardize, _validate
from .query_manager import _lookup, _search
from .sqlrecord import (
    Registry,
    SQLRecord,
    _adjust_is_latest_when_deleting_is_versioned,
)

if TYPE_CHECKING:
    from bionty.models import (
        CellLine,
        CellMarker,
        CellType,
        DevelopmentalStage,
        Disease,
        Ethnicity,
        ExperimentalFactor,
        Gene,
        Organism,
        Pathway,
        Phenotype,
        Protein,
        Tissue,
    )
    from pertdb.models import (
        Biologic,
        CombinationPerturbation,
        Compound,
        CompoundPerturbation,
        EnvironmentalPerturbation,
        GeneticPerturbation,
        PerturbationTarget,
    )

    from lamindb.base.types import ListLike, StrField
    from lamindb.models import (
        Artifact,
        Branch,
        Collection,
        Feature,
        Project,
        Record,
        Reference,
        Run,
        Schema,
        Space,
        Storage,
        Transform,
        ULabel,
        User,
    )

T = TypeVar("T")


pd.set_option("display.max_columns", 200)


def get_keys_from_df(data: list, registry: SQLRecord) -> list[str]:
    if len(data) > 0:
        if isinstance(data[0], dict):
            keys = list(data[0].keys())
        else:
            keys = list(data[0].__dict__.keys())
            if "_state" in keys:
                keys.remove("_state")
    else:
        keys = [
            field.name
            for field in registry._meta.fields
            if not isinstance(field, models.ForeignKey)
        ]
        keys += [
            f"{field.name}_id"
            for field in registry._meta.fields
            if isinstance(field, models.ForeignKey)
        ]
    return keys


def get_default_branch_ids(branch: Branch | None = None) -> list[int]:
    """Return branch IDs to include in default queries.

    By default, queries include records on the main branch (branch_id=1) but exclude trashed (branch_id=-1)
    and archived records (branch_id=0). This matches behavior of familiar tools like GitHub, Slack, and
    email clients.

    If a user switches to another branch via `lamin switch branch`, the main branch will still be included.

    Returns:
        List containing the default branch and current branch if different.
    """
    if branch is None:
        branch_id = setup_settings.branch.id
    else:
        branch_id = branch.id
    branch_ids = [branch_id]
    if branch_id != 1:  # add the main branch by default
        branch_ids.append(1)
    return branch_ids


def one_helper(
    self: QuerySet | SQLRecordList,
    does_not_exist_msg: str | None = None,
    raise_doesnotexist: bool = True,
    not_exists: bool | None = None,
    raise_multipleresultsfound: bool = True,
):
    if not_exists is None:
        if isinstance(self, SQLRecordList):
            not_exists = len(self) == 0
        else:
            not_exists = not self.exists()  # type: ignore
    if not_exists:
        if raise_doesnotexist:
            raise DoesNotExist(does_not_exist_msg)
        else:
            return None
    elif len(self) > 1:
        if raise_multipleresultsfound:
            raise MultipleResultsFound(self)
        else:
            return self[0]
    else:
        return self[0]


def get_backward_compat_filter_kwargs(queryset, expressions):
    from lamindb.models import (
        Artifact,
        Feature,
    )

    if issubclass(queryset.model, IsVersioned):
        name_mappings = {
            "version": "version_tag",
        }
    else:
        name_mappings = {}

    if queryset.model is Artifact:
        name_mappings.update(
            {
                "transform": "run__transform",
                "feature_sets": "schemas",
            }
        )
    if queryset.model is Feature:
        name_mappings.update(
            {
                "dtype": "_dtype_str",
                "dtype_as_str": "_dtype_str",
            }
        )

    # If no mappings to apply, return expressions as-is
    if not name_mappings:
        return expressions
    was_list = False
    if isinstance(expressions, list):
        was_list = True
        expressions = {field: True for field in expressions}
    mapped = {}
    for field, value in expressions.items():
        parts = field.split("__")
        if parts[0] in name_mappings:
            # Issue deprecation warnings
            if queryset.model is Artifact and parts[0] == "feature_sets":
                warnings.warn(
                    "Querying Artifact by `feature_sets` is deprecated. Use `schemas` instead.",
                    DeprecationWarning,
                    stacklevel=4,
                )
            elif queryset.model is Feature and parts[0] == "dtype":
                warnings.warn(
                    "Querying Feature by `dtype` is deprecated. Use `dtype_as_str` instead. "
                    "Notice the new dtype encoding format for Record and ULabel subtypes.",
                    DeprecationWarning,
                    stacklevel=4,
                )
            new_field = name_mappings[parts[0]] + (
                "__" + "__".join(parts[1:]) if len(parts) > 1 else ""
            )
            mapped[new_field] = value
        else:
            mapped[field] = value
    return list(mapped.keys()) if was_list else mapped


def process_expressions(queryset: QuerySet, queries: tuple, expressions: dict) -> dict:
    def _map_databases(value: Any, key: str, target_db: str) -> tuple[str, Any]:
        if isinstance(value, SQLRecord):
            if value._state.db != target_db:
                logger.warning(
                    f"passing record from database {value._state.db} to query {target_db}, matching on uid '{value.uid}'"
                )
                return f"{key}__uid", value.uid
            return key, value

        if (
            key.endswith("__in")
            and isinstance(value, IterableType)
            and not isinstance(value, str)
        ):
            if any(
                isinstance(v, SQLRecord) and v._state.db != target_db for v in value
            ):
                logger.warning(
                    f"passing records from another database to query {target_db}, matching on uids"
                )
                return key.replace("__in", "__uid__in"), [
                    v.uid if isinstance(v, SQLRecord) else v for v in value
                ]
            return key, value

        return key, value

    branch_fields = {"branch", "branch_id"}
    branch_prefixes = ("branch__", "branch_id__")

    def queries_contain_branch(queries: tuple) -> bool:
        """Check if any Q object in queries references branch or branch_id."""

        def check_q_object(q: Q) -> bool:
            # Q objects store their conditions in q.children
            for child in q.children:
                if isinstance(child, tuple) and len(child) == 2:
                    # Normal condition: (key, value)
                    key = child[0]
                    if key in branch_fields or key.startswith(branch_prefixes):
                        return True
                elif isinstance(child, Q):
                    # Nested Q object
                    if check_q_object(child):
                        return True
            return False

        return any(check_q_object(q) for q in queries if isinstance(q, Q))

    expressions = get_backward_compat_filter_kwargs(
        queryset,
        expressions,
    )
    if issubclass(queryset.model, SQLRecord):
        # branch_id is set to 1 unless expressions contains id, uid or hash
        id_uid_hash = {"id", "uid", "hash", "id__in", "uid__in", "hash__in"}
        if not any(expression in id_uid_hash for expression in expressions):
            expressions_have_branch = False
            for expression in expressions:
                if expression in branch_fields or expression.startswith(
                    branch_prefixes
                ):
                    expressions_have_branch = True
                    break
            if not expressions_have_branch and not queries_contain_branch(queries):
                expressions["branch_id__in"] = get_default_branch_ids()
            else:
                # if branch_id is None, do not apply a filter
                # otherwise, it would mean filtering for NULL values, which doesn't make
                # sense for a non-NULLABLE column
                if "branch_id" in expressions and expressions["branch_id"] is None:
                    expressions.pop("branch_id")
                if "branch" in expressions and expressions["branch"] is None:
                    expressions.pop("branch")

    if queryset._db is not None:
        # only check for database mismatch if there is a defined database on the
        # queryset
        return dict(
            (
                _map_databases(value, key, queryset._db)
                for key, value in expressions.items()
            )
        )
    else:
        return expressions


def get(
    registry_or_queryset: Registry | BasicQuerySet,
    idlike: int | str | None = None,
    **expressions,
) -> SQLRecord:
    if isinstance(registry_or_queryset, BasicQuerySet):
        # not QuerySet but only BasicQuerySet
        assert not isinstance(registry_or_queryset, QuerySet)  # noqa: S101

        qs = registry_or_queryset
        registry = qs.model
    else:
        qs = BasicQuerySet(model=registry_or_queryset)
        registry = registry_or_queryset

    if isinstance(idlike, int):
        return qs.get(id=idlike)
    elif isinstance(idlike, str):
        NAME_FIELD = (
            registry._name_field if hasattr(registry, "_name_field") else "name"
        )
        DOESNOTEXIST_MSG = f"No record found with uid '{idlike}'. Did you forget a keyword as in {registry.__name__}.get({NAME_FIELD}='{idlike}')?"
        # this is the case in which the user passes an under-specified uid
        if issubclass(registry, IsVersioned) and len(idlike) <= registry._len_stem_uid:
            new_qs = qs.filter(uid__startswith=idlike, is_latest=True)
            not_exists = None
            if not new_qs.exists():
                # also try is_latest is False due to nothing found
                new_qs = qs.filter(uid__startswith=idlike, is_latest=False)
            else:
                not_exists = False
            # it doesn't make sense to raise MultipleResultsFound when querying with an
            # underspecified uid
            return one_helper(
                new_qs,
                DOESNOTEXIST_MSG,
                not_exists=not_exists,
                raise_multipleresultsfound=False,
            )
        else:
            qs = qs.filter(uid__startswith=idlike)
            return one_helper(qs, DOESNOTEXIST_MSG)
    else:
        assert idlike is None  # noqa: S101
        expressions = process_expressions(qs, [], expressions)
        # inject is_latest for consistency with idlike
        is_latest_was_not_in_expressions = "is_latest" not in expressions
        if issubclass(registry, IsVersioned) and is_latest_was_not_in_expressions:
            expressions["is_latest"] = True
        try:
            return qs.get(**expressions)
        except registry.DoesNotExist as e:
            # handle the case in which the is_latest injection led to a missed query
            if "is_latest" in expressions and is_latest_was_not_in_expressions:
                expressions.pop("is_latest")
                result = qs.filter(**expressions).order_by("-created_at").first()
                if result is not None:
                    return result
            raise e


class SQLRecordList(UserList, Generic[T]):
    """Is ordered, can't be queried, but has `.to_dataframe()`."""

    def __init__(self, records: Iterable[T]):
        if isinstance(records, list):
            self.data = records  # Direct assignment if already a list, no copy
        else:
            super().__init__(records)  # Let UserList handle the conversion

    def to_dataframe(self) -> pd.DataFrame:
        keys = get_keys_from_df(self.data, self.data[0].__class__)
        values = [record.__dict__ for record in self.data]
        return pd.DataFrame(values, columns=keys)

    @deprecated(new_name="to_dataframe")
    def df(self) -> pd.DataFrame:
        return self.to_dataframe()

    def to_list(
        self, field: str | None = None
    ) -> list[str]:  # meaningful to be parallel with to_list() in QuerySet
        if field is None:
            return self.data
        return [getattr(record, field) for record in self.data]

    def one(self) -> T:
        """Exactly one result. Throws error if there are more or none."""
        return one_helper(self)

    def save(self) -> SQLRecordList[T]:
        """Save all records to the database."""
        from lamindb.models.save import save

        save(self)
        return self


def get_basic_field_names(
    qs: QuerySet,
    include: list[str],
    features_input: bool | list[str] | str,
) -> list[str]:
    exclude_field_names = ["updated_at"]
    include_private_fields = False
    if "privates" in include:
        include_private_fields = True
        include.remove("privates")
    field_names = [
        field.name
        for field in qs.model._meta.fields
        if (
            not isinstance(field, models.ForeignKey)
            and field.name not in exclude_field_names
            and (
                not field.name.startswith("_")
                or include_private_fields
                or (field.name == "_dtype_str" and qs.model.__name__ == "Feature")
            )
        )
    ]
    for field_name in [
        "version_tag",
        "is_latest",
        "is_locked",
        "is_type",
        "created_at",
        "updated_at",
    ]:
        if field_name in field_names:
            field_names.append(field_names.pop(field_names.index(field_name)))
    field_names += [
        f"{field.name}_id"
        for field in qs.model._meta.fields
        if isinstance(field, models.ForeignKey)
    ]
    # move uid to first position if present
    if "uid" in field_names:
        field_names.insert(0, field_names.pop(field_names.index("uid")))

    # move primary key to second position if present
    pk = qs.model._meta.pk.name if qs.model._meta.pk else None
    if pk and pk in field_names:
        field_names.insert(1, field_names.pop(field_names.index(pk)))
    if (
        include or features_input
    ):  # if there is features_input, reduce fields to just the first 3
        subset_field_names = field_names[:3]
        intersection = set(field_names) & set(include)
        subset_field_names += list(intersection)
        field_names = subset_field_names
    return field_names


def get_feature_annotate_kwargs(
    registry: Registry,
    features: bool | list[str] | str | None,
    qs: QuerySet | None = None,
) -> tuple[dict[str, Any], QuerySet, dict[str, Any]]:
    from lamindb.models import (
        Artifact,
        Feature,
        Record,
        RecordJson,
        ULabel,
    )
    from lamindb.models.feature import parse_dtype

    if registry not in {Artifact, Record}:
        raise ValueError(
            f"features=True is only applicable for Artifact and Record, not {registry.__name__}"
        )

    feature_ids = []
    if features == "queryset":
        ids_list = qs.values_list("id", flat=True)
        for obj in registry._meta.related_objects:
            related_name_attr = getattr(registry, obj.related_name, None)
            if related_name_attr is None or not hasattr(related_name_attr, "through"):
                continue
            link_model = related_name_attr.through
            if (
                not hasattr(link_model, "feature")
                or link_model.__name__ == "Record_parents"
            ):
                continue
            filter_field = registry.__name__.lower()
            if not hasattr(link_model, filter_field):
                potential_fields = []
                for field in link_model._meta.get_fields():
                    if field.is_relation and field.related_model is registry:
                        potential_fields.append(field.name)
                if len(potential_fields) == 1:
                    filter_field = potential_fields[0]
                else:
                    continue
            links = link_model.objects.using(qs.db).filter(
                **{filter_field + "_id__in": ids_list}
            )
            feature_ids_for_link_model = links.values_list("feature__id", flat=True)
            feature_ids += feature_ids_for_link_model
        if registry is Record:
            # this request is not strictly necessary, but it makes the resulting reshaped
            # dataframe consistent
            feature_ids += RecordJson.filter(record_id__in=ids_list).values_list(
                "feature__id", flat=True
            )
        feature_ids = list(set(feature_ids))  # remove duplicates

    feature_qs = Feature.connect(None if qs is None else qs.db).filter(
        _dtype_str__isnull=False
    )
    if isinstance(features, list):
        feature_qs = feature_qs.filter(name__in=features)
        if len(features) != feature_qs.count():
            logger.warning(
                f"found features and passed features differ:\n - passed: {features}\n - found: {feature_qs.to_list('name')}"
            )
    elif feature_ids:
        feature_qs = feature_qs.filter(id__in=feature_ids)
    else:
        feature_qs = feature_qs.filter(
            ~Q(_dtype_str__startswith="cat[")
            | Q(_dtype_str__startswith="cat[ULabel")
            | Q(_dtype_str__startswith="cat[Record")
        )
        logger.important(
            f"queried for all categorical features of dtypes Record or ULabel and non-categorical features: ({len(feature_qs)}) {feature_qs.to_list('name')}"
        )
    # Get the categorical features
    cat_feature_types = {
        parse_dtype(feature._dtype_str)[0]["registry_str"]
        for feature in feature_qs
        if feature._dtype_str.startswith("cat[")
        or feature._dtype_str.startswith("list[cat[")
    }
    # fields to annotate
    cat_feature_fields = defaultdict(list)
    for feature in feature_qs:
        dtype_str = feature._dtype_str
        if dtype_str.startswith("cat[") or dtype_str.startswith("list[cat["):
            dtype_info = parse_dtype(dtype_str)[0]
            registry_str = dtype_info["registry_str"]
            field_name = dtype_info["field_str"]
            cat_feature_fields[registry_str].append(field_name)
    # Get relationships of labels and features
    link_models_on_models = {
        getattr(
            registry, obj.related_name
        ).through.__get_name_with_module__(): obj.related_model
        for obj in registry._meta.related_objects
        if obj.related_model.__get_name_with_module__() in cat_feature_types
        and hasattr(getattr(registry, obj.related_name), "through")
        and hasattr(getattr(registry, obj.related_name).through, "feature_id")
    }
    if registry is Artifact:
        link_models_on_models["ArtifactULabel"] = ULabel
    else:
        link_models_on_models["RecordRecord"] = Record
    link_attributes_on_models = {
        obj.related_name: link_models_on_models[
            obj.related_model.__get_name_with_module__()
        ]
        for obj in registry._meta.related_objects
        if (
            obj.related_model.__get_name_with_module__() in link_models_on_models
            and (
                not obj.related_name.startswith("links_record")
                if registry is Record
                else True
            )
        )
    }
    # Prepare Django's annotate for features with filtering
    filtered_relations = {}
    annotate_kwargs = {}

    for link_attr, feature_type_model in link_attributes_on_models.items():
        feature_type = feature_type_model.__get_name_with_module__()
        if link_attr == "links_project" and registry is Record:
            # we're only interested in _values_project when "annotating" records
            continue

        # Determine field name
        if registry is Artifact:
            field_name = (
                feature_type.split(".")[1] if "." in feature_type else feature_type
            ).lower()
        else:
            field_name = "value"

        # Determine if this value model needs branch filtering
        # Skip user relations (RecordUser, ArtifactUser don't have branch)
        should_filter_branch = link_attr not in {"values_user", "links_user"}

        # Create filtered relation for the value model
        value_relation_path = f"{link_attr}__{field_name}"
        filtered_value_relation_name = f"filtered_{link_attr}_{field_name}"

        if should_filter_branch:
            filtered_relations[filtered_value_relation_name] = FilteredRelation(
                value_relation_path,
                condition=Q(
                    **{
                        f"{value_relation_path}__branch_id__in": get_default_branch_ids()
                    }
                ),
            )
        else:
            # No branch filtering needed
            filtered_relations[filtered_value_relation_name] = FilteredRelation(
                value_relation_path
            )

        # Add annotation for feature name (feature doesn't have branch_id)
        annotate_kwargs[f"{link_attr}__feature__name"] = F(
            f"{link_attr}__feature__name"
        )

        # Add annotations for categorical feature fields using the filtered relation
        for field in cat_feature_fields[feature_type]:
            annotate_kwargs[f"{link_attr}__{field_name}__{field}"] = F(
                f"{filtered_value_relation_name}__{field}"
            )

    # Handle JSON values (no branch filtering needed)
    json_values_attribute = "json_values" if registry is Artifact else "values_json"
    annotate_kwargs[f"{json_values_attribute}__feature__name"] = F(
        f"{json_values_attribute}__feature__name"
    )
    annotate_kwargs[f"{json_values_attribute}__value"] = F(
        f"{json_values_attribute}__value"
    )

    return annotate_kwargs, feature_qs, filtered_relations


# https://claude.ai/share/16280046-6ae5-4f6a-99ac-dec01813dc3c
def analyze_lookup_cardinality(
    model_class: SQLRecord, lookup_paths: list[str] | None
) -> dict[str, str]:
    """Analyze lookup cardinality.

    Analyzes Django model lookups to determine if they will result in
    one-to-one or one-to-many relationships when used in annotations.

    Args:
        model_class: The Django model class to analyze
        include: List of lookup paths (e.g. ["created_by__name", "ulabels__name"])

    Returns:
        Dictionary mapping lookup paths to either 'one' or 'many'
    """
    result = {}  # type: ignore
    if lookup_paths is None:
        return result
    for lookup_path in lookup_paths:
        parts = lookup_path.split("__")
        current_model = model_class
        is_many = False

        # Walk through each part of the lookup path
        for part in parts[:-1]:  # Exclude the last part as it's an attribute
            field = None

            # Handle reverse relations
            for f in current_model._meta.get_fields():
                if isinstance(f, ForeignObjectRel) and f.get_accessor_name() == part:
                    field = f
                    is_many = not f.one_to_one
                    if hasattr(f, "field"):
                        current_model = f.field.model
                    break

            # Handle forward relations
            if field is None:
                field = current_model._meta.get_field(part)
                if isinstance(field, ManyToManyField):
                    is_many = True
                    current_model = field.remote_field.model
                elif isinstance(field, ForeignKey):
                    current_model = field.remote_field.model

        result[lookup_path] = "many" if is_many else "one"

    return result


def reorder_subset_columns_in_df(
    df: pd.DataFrame, column_order: list[str], position=3
) -> pd.DataFrame:
    """Reorder subset of columns in dataframe to specified position."""
    valid_columns = [col for col in column_order if col in df.columns]
    all_cols = df.columns.tolist()
    remaining_cols = [col for col in all_cols if col not in valid_columns]
    new_order = remaining_cols[:position] + valid_columns + remaining_cols[position:]
    return df[new_order]


def encode_lamindb_fields_as_columns(
    registry: Registry, fields: str | list[str]
) -> str | dict[str, str]:
    """Encode laminDB specific fields in dataframe with __lamindb_{model_name}_{field_name}__.

    This is needed when reshaping dataframes with features to avoid conflicts between
    laminDB fields and feature names.
    """

    def encode(field: str) -> str:
        return f"__lamindb_{registry._meta.model_name}_{field}__"

    registry_field_names = {field.name for field in registry._meta.concrete_fields}

    if isinstance(fields, str):
        return encode(fields) if fields in registry_field_names else fields

    return {field: encode(field) for field in fields if field in registry_field_names}


# https://lamin.ai/laminlabs/lamindata/transform/BblTiuKxsb2g0003
# https://claude.ai/chat/6ea2498c-944d-4e7a-af08-29e5ddf637d2
def reshape_annotate_result(
    registry: Registry,
    df: pd.DataFrame,
    field_names: list[str],
    cols_from_include: dict[str, str] | None,
    feature_qs: QuerySet | None,
) -> pd.DataFrame:
    """Reshapes tidy table to wide format.

    Args:
        registry: The registry model (e.g., Artifact)
        df: Input dataframe with experimental data
        field_names: List of basic fields to include in result
        cols_from_include: Dict specifying additional columns to process with types
            ('one' or 'many'), e.g., {'ulabels__name': 'many', 'created_by__name': 'one'}
        feature_qs: QuerySet of features
    """
    from lamindb.models import Artifact

    cols_from_include = cols_from_include or {}

    # Initialize result with basic fields (need a copy since we're modifying it)
    result = df[field_names].copy()
    pk_name = registry._meta.pk.name

    # ========== no features requested ==========
    if feature_qs is None or not feature_qs.exists():
        if cols_from_include:
            result = process_cols_from_include(df, result, cols_from_include, pk_name)
        return result.drop_duplicates(subset=[pk_name])

    # ========== process features ==========

    # Encode Django field names to avoid conflicts with feature names
    fields_map = encode_lamindb_fields_as_columns(registry, df.columns)
    df_encoded = df.rename(columns=fields_map)
    result_encoded = result.rename(columns=fields_map)
    pk_name_encoded = fields_map.get(pk_name)  # type: ignore

    # --- Process JSON-stored feature values ---
    json_values_attribute = "json_values" if registry is Artifact else "values_json"
    feature_name_col = f"{json_values_attribute}__feature__name"
    feature_value_col = f"{json_values_attribute}__value"

    if all(col in df_encoded.columns for col in [feature_name_col, feature_value_col]):
        # Separate dict and non-dict values for different aggregation strategies
        is_dict_or_list = df_encoded[feature_value_col].apply(
            lambda x: isinstance(x, (dict, list))
        )
        dict_or_list_df = df_encoded[is_dict_or_list]
        non_dict_or_list_df = df_encoded[~is_dict_or_list]

        # Aggregate: sets for non-dict values, first for dict values
        groupby_cols = [pk_name_encoded, feature_name_col]
        non_dict_or_list_features = non_dict_or_list_df.groupby(groupby_cols)[
            feature_value_col
        ].agg(set)
        dict_or_list_features = dict_or_list_df.groupby(groupby_cols)[
            feature_value_col
        ].agg("first")

        # Combine and pivot to wide format
        combined_features = pd.concat(
            [non_dict_or_list_features, dict_or_list_features]
        )
        feature_values = combined_features.unstack().reset_index()

        if not feature_values.empty:
            result_encoded = result_encoded.join(
                feature_values.set_index(pk_name_encoded),
                on=pk_name_encoded,
            )

    # --- Process categorical/linked features ---
    links_prefix = "links_" if registry is Artifact else ("links_", "values_")
    links_features = [
        col
        for col in df.columns
        if "feature__name" in col and col.startswith(links_prefix)
    ]

    if links_features:
        result_encoded = process_links_features(
            df_encoded,
            result_encoded,
            links_features,
            feature_qs,
            pk_name_encoded,
        )

    # --- Apply type conversions based on feature metadata ---
    def extract_and_check_scalar(series: pd.Series) -> tuple[pd.Series, bool]:
        """Extract single elements and return if column is now scalar."""
        has_multiple_values = False

        def extract_and_track(value):
            nonlocal has_multiple_values
            if not hasattr(value, "__len__") or isinstance(value, str):
                return value
            if len(value) != 1:
                has_multiple_values = True
                return value
            return next(iter(value))

        extracted = series.apply(extract_and_track)
        is_scalar = not has_multiple_values
        return extracted, is_scalar

    for feature in feature_qs:
        if feature.name not in result_encoded.columns:
            continue

        result_encoded[feature.name], is_scalar = extract_and_check_scalar(
            result_encoded[feature.name]
        )

        if is_scalar:
            dtype_str = feature._dtype_str
            if dtype_str.startswith("cat"):
                result_encoded[feature.name] = result_encoded[feature.name].astype(
                    "category"
                )
            if dtype_str == "datetime":
                # format and utc args are needed for mixed data
                # pandera expects timezone-naive datetime objects, and hence,
                # we need to localize with None
                result_encoded[feature.name] = pd.to_datetime(
                    result_encoded[feature.name], format="ISO8601", utc=True
                ).dt.tz_localize(None)
            if dtype_str == "date":
                # see comments for datetime
                result_encoded[feature.name] = (
                    pd.to_datetime(
                        result_encoded[feature.name],
                        format="ISO8601",
                        utc=True,
                    )
                    .dt.tz_localize(None)
                    .dt.date
                )
            if dtype_str == "bool":
                result_encoded[feature.name] = result_encoded[feature.name].astype(
                    "boolean"
                )

        dtype_str = feature._dtype_str
        if dtype_str.startswith("list"):
            mask = result_encoded[feature.name].notna()
            result_encoded.loc[mask, feature.name] = result_encoded.loc[
                mask, feature.name
            ].apply(lambda x: list(x) if isinstance(x, (set, list)) else [x])

        if dtype_str == "dict":
            # this is the case when a dict is stored as a string; won't happen
            # within lamindb but might for external data
            if isinstance(result_encoded[feature.name].iloc[0], str):
                result_encoded[feature.name] = result_encoded[feature.name].apply(
                    lambda x: ast.literal_eval(x) if isinstance(x, str) else x
                )

    # --- Finalize result ---

    # Reorder columns to prioritize features
    result_encoded = reorder_subset_columns_in_df(
        result_encoded,
        feature_qs.to_list("name"),  # type: ignore
    )

    # Process additional included columns
    if cols_from_include:
        cols_from_include_encoded = {
            fields_map.get(k, k): v  # type: ignore
            for k, v in cols_from_include.items()
        }
        result_encoded = process_cols_from_include(
            df_encoded, result_encoded, cols_from_include_encoded, pk_name_encoded
        )

    # Decode field names back to original, except where conflicts exist
    # (e.g., if a feature is also named 'id', keep the encoded field name)
    decode_map = {
        encoded: original
        for original, encoded in fields_map.items()  # type: ignore
        if original not in result_encoded.columns
    }

    return result_encoded.drop_duplicates(subset=[pk_name_encoded]).rename(
        columns=decode_map
    )


def process_links_features(
    df: pd.DataFrame,
    result: pd.DataFrame,
    feature_cols: list[str],
    feature_qs: QuerySet | None,
    pk_name: str = "id",
) -> pd.DataFrame:
    """Process links_XXX feature columns."""
    from lamindb.models.feature import parse_dtype

    # this loops over different entities that might be linked under a feature
    for feature_col in feature_cols:
        links_attribute = "links_" if feature_col.startswith("links_") else "values_"
        regex = f"{links_attribute}(.+?)__feature__name"
        prefix = re.match(regex, feature_col).group(1)

        value_cols = [
            col
            for col in df.columns
            if col.startswith(f"{links_attribute}{prefix}__")
            and "feature__name" not in col
        ]

        if not value_cols:
            continue

        value_col = value_cols[0]
        feature_names = df[feature_col].unique()
        feature_names = feature_names[~pd.isna(feature_names)]

        for feature in feature_qs:
            if feature.name not in feature_names:
                continue
            if feature.name in result.columns:
                continue
            field_name = parse_dtype(feature._dtype_str)[0]["field_str"]
            value_col = [c for c in value_cols if c.endswith(f"__{field_name}")][0]
            mask = (df[feature_col] == feature.name) & df[value_col].notna()
            feature_values = df[mask].groupby(pk_name)[value_col].agg(set)
            result.insert(3, feature.name, result[pk_name].map(feature_values))

    return result


def process_cols_from_include(
    df: pd.DataFrame,
    result: pd.DataFrame,
    extra_columns: dict[str, str],
    pk_name: str = "id",
) -> pd.DataFrame:
    """Process additional columns based on their specified types."""
    for col, col_type in extra_columns.items():
        if col not in df.columns:
            continue
        if col in result.columns:
            continue

        values = df.groupby(pk_name)[col].agg(set if col_type == "many" else "first")
        result.insert(3, col, result[pk_name].map(values))

    return result


def _queryset_class_factory(
    registry: Registry, queryset_cls: type[models.QuerySet]
) -> type[models.QuerySet]:
    from lamindb.models import Artifact, ArtifactSet

    # If the model is Artifact, create a new class for BasicQuerySet or QuerySet that inherits from ArtifactSet.
    # This allows to add artifact specific functionality to all classes inheriting from BasicQuerySet.
    # Thus all query sets of artifacts (and only of artifacts) will have functions from ArtifactSet.
    if registry is Artifact and not issubclass(queryset_cls, ArtifactSet):
        new_cls = type(
            "Artifact" + queryset_cls.__name__, (queryset_cls, ArtifactSet), {}
        )
    else:
        new_cls = queryset_cls
    return new_cls


class BasicQuerySet(models.QuerySet):
    """Sets of records returned by queries.

    See Also:

        `django QuerySet <https://docs.djangoproject.com/en/stable/ref/models/querysets/>`__

    Examples:

        Any filter statement produces a query set::

            queryset = Registry.filter(name__startswith="keyword")
    """

    def __new__(cls, model=None, query=None, using=None, hints=None):
        # see comments in _queryset_class_factory
        return object.__new__(_queryset_class_factory(model, cls))

    def _to_class(
        self, cls: type[models.QuerySet], copy: bool = True
    ) -> models.QuerySet:
        qs = self.all() if copy else self
        qs.__class__ = cls
        return qs

    def _to_basic(self, copy: bool = True) -> BasicQuerySet:
        cls = _queryset_class_factory(self.model, BasicQuerySet)
        return self._to_class(cls, copy)

    def _to_non_basic(self, copy: bool = True) -> QuerySet:
        cls = _queryset_class_factory(self.model, QuerySet)
        return self._to_class(cls, copy)

    @doc_args(SQLRecord.to_dataframe.__doc__)
    def to_dataframe(
        self,
        *,
        include: str | list[str] | None = None,
        features: str | list[str] | None = None,
        limit: int | None = 100,
        order_by: str | None = "-id",
    ) -> pd.DataFrame:
        """{}"""  # noqa: D415
        if (
            self.model.__name__ == "Artifact"
            and "kind" not in str(self.query.where)
            and self.query.low_mark
            == 0  # this should be 0, not None, it represent OFFSET = 0
            and self.query.high_mark
            is None  # this should be None, it represent _no_ LIMIT
        ):
            subset = self.exclude(**{"kind__startswith": "__lamindb"})
        else:
            subset = self
        # check if queryset is already ordered
        is_ordered = bool(subset.query.order_by)
        # Only apply order_by if not already ordered and order_by is specified
        if not is_ordered and order_by is not None:
            subset = subset.order_by(order_by)
        if limit is not None:
            subset = subset[:limit]
        if include is None:
            include_input = []
        elif isinstance(include, str):
            include_input = [include]
        else:
            include_input = include
        if "features" in include_input:
            include_input.remove("features")
            if features is None:
                # indicate the default features with True
                # should refactor this in the future
                features = True  # type: ignore
        features_input = [] if features is None else features
        include = get_backward_compat_filter_kwargs(subset, include_input)
        field_names = get_basic_field_names(subset, include_input, features_input)

        annotate_kwargs = {}
        filtered_relations = {}  # type: ignore
        feature_qs = None
        if features:
            feature_annotate_kwargs, feature_qs, filtered_relations = (
                get_feature_annotate_kwargs(subset.model, features, subset)
            )
            annotate_kwargs.update(feature_annotate_kwargs)
        if include_input:
            include_input = include_input.copy()[::-1]  # type: ignore
            include_kwargs = {s: F(s) for s in include_input if s not in field_names}
            annotate_kwargs.update(include_kwargs)
        if annotate_kwargs:
            id_subquery = subset.values("id")
            # for annotate, we want the queryset without filters so that joins don't affect the annotations
            query_set_without_filters = subset.model.objects.using(subset.db).filter(
                id__in=Subquery(id_subquery)
            )
            if subset.query.order_by:
                # Apply the same ordering to the new queryset
                query_set_without_filters = query_set_without_filters.order_by(
                    *subset.query.order_by
                )
            if filtered_relations:
                query_set_without_filters = query_set_without_filters.annotate(
                    **filtered_relations
                )
            queryset = query_set_without_filters.annotate(**annotate_kwargs)
        else:
            queryset = subset

        # our main problem with this approach is that we lose ordering in categorical lists
        # we'd need to respect ordering through the primary key on the link table, but that's
        # another refactoring effort
        # we have the correct ordering in `features.get_values()`, though
        df = pd.DataFrame(queryset.values(*field_names, *list(annotate_kwargs.keys())))
        if len(df) == 0:
            df = pd.DataFrame({}, columns=field_names)
            return df
        cols_from_include = analyze_lookup_cardinality(self.model, include_input)  # type: ignore
        df_reshaped = reshape_annotate_result(
            self.model, df, field_names, cols_from_include, feature_qs
        )
        pk_name = self.model._meta.pk.name
        encoded_pk_name = encode_lamindb_fields_as_columns(self.model, pk_name)
        if encoded_pk_name in df_reshaped.columns:
            df_reshaped = df_reshaped.set_index(encoded_pk_name)
        else:
            pk_column_name = pk_name if pk_name in df.columns else f"{pk_name}_id"
            if pk_column_name in df_reshaped.columns:
                df_reshaped = df_reshaped.set_index(pk_column_name)

        # cast floats and ints where appropriate
        # this is currently needed because the UI writes into the JSON field through JS
        # and thus a `10` might be a float, not an int
        # note: also type casting within reshape_annotate_result
        if feature_qs is not None:
            for feature in feature_qs:
                if feature.name in df_reshaped.columns:
                    current_dtype = df_reshaped[feature.name].dtype
                    dtype_str = feature._dtype_str
                    if dtype_str == "int" and not pd.api.types.is_integer_dtype(
                        current_dtype
                    ):
                        df_reshaped[feature.name] = df_reshaped[feature.name].astype(
                            "Int64"  # nullable integer dtype
                        )
                    elif dtype_str == "float" and not pd.api.types.is_float_dtype(
                        current_dtype
                    ):
                        df_reshaped[feature.name] = df_reshaped[feature.name].astype(
                            float
                        )
        return df_reshaped

    @deprecated(new_name="to_dataframe")
    def df(
        self,
        include: str | list[str] | None = None,
        features: bool | list[str] | str | None = None,
    ) -> pd.DataFrame:
        return self.to_dataframe(include=include, features=features)

    def describe(self, return_str: bool = False) -> str | None:
        """Describe the query set to learn about available fields."""
        return self.model.describe(return_str=return_str)

    def delete(self, *args, permanent: bool | None = None, **kwargs):
        """Delete all records in the query set.

        Args:
            permanent: Whether to permanently delete the record (skips trash).
                Is only relevant for records that have the `branch` field.
                If `None`, uses soft delete for records that have the `branch` field, hard delete otherwise.

        Note:
            Calling `delete()` twice on the same queryset does NOT permanently delete in bulk operations.
            Use `permanent=True` for actual deletion.

        Examples:

            For a `QuerySet` object `qs`, call::

                qs.delete()
        """
        from lamindb.models import Artifact, Collection, Run, Storage, Transform

        if self.model is Run:
            if permanent is True:
                from .run import _permanent_delete_runs

                _permanent_delete_runs(self)
                return
            if permanent is not True:
                self.update(branch_id=-1)
                return
        if self.model is Transform:
            if permanent is True:
                from .transform import _permanent_delete_transforms

                _permanent_delete_transforms(self)
                return
            if permanent is not True:
                transform_ids = list(self.values_list("pk", flat=True))
                if transform_ids:
                    _adjust_is_latest_when_deleting_is_versioned(
                        Transform, self.db or "default", transform_ids
                    )
                self.update(branch_id=-1, is_latest=False)
                return
        # Artifact, Collection: non-trivial delete behavior, handle in a loop
        if self.model in {Artifact, Collection}:
            for record in self:
                record.delete(*args, permanent=permanent, **kwargs)
        elif self.model is Storage:  # storage does not have soft delete
            if permanent is False:
                raise ValueError(
                    "Soft delete is not possible for Storage, "
                    "use 'permanent=True' or 'permanent=None' for permanent deletion."
                )
            for record in self:
                record.delete()
        else:
            if not permanent and hasattr(self.model, "branch_id"):
                logger.warning("moved records to trash (branch_id = -1)")
                self.update(branch_id=-1)
            else:
                if permanent is False:
                    raise ValueError(
                        f"Soft delete is not possible for {self.model.__name__}, "
                        "use 'permanent=True' for permanent deletion."
                    )
                super().delete(*args, **kwargs)

    def to_list(self, field: str | None = None) -> list[SQLRecord] | list[str]:
        """Populate an (unordered) list with the results.

        Note that the order in this list is only meaningful if you ordered the underlying query set with `.order_by()`.

        Examples::

            queryset.to_list()  # list of records
            queryset.to_list("name")  # list of values
        """
        if field is None:
            return list(self)
        else:
            # list casting is necessary because values_list does not return a list
            return list(self.values_list(field, flat=True))

    def first(self) -> SQLRecord | None:
        """If non-empty, the first result in the query set, otherwise ``None``.

        Examples::

            queryset.first()
        """
        if len(self) == 0:
            return None
        return self[0]

    def one(self) -> SQLRecord:
        """Exactly one result. Raises error if there are more or none."""
        return one_helper(self)

    def one_or_none(self) -> SQLRecord | None:
        """At most one result. Returns it if there is one, otherwise returns ``None``.

        Examples::

            ULabel.filter(name="benchmark").one_or_none()
            ULabel.filter(name="non existing label").one_or_none()
        """
        return one_helper(self, raise_doesnotexist=False)

    @doc_args(_search.__doc__)
    def search(self, string: str, **kwargs):
        """{}"""  # noqa: D415
        return _search(cls=self, string=string, **kwargs)

    @doc_args(_lookup.__doc__)
    def lookup(self, field: StrField | None = None, **kwargs) -> NamedTuple:
        """{}"""  # noqa: D415
        return _lookup(cls=self, field=field, **kwargs)

    # -------------------------------------------------------------------------------------
    # CanCurate
    # -------------------------------------------------------------------------------------

    @doc_args(CanCurate.validate.__doc__)
    def validate(self, values: ListLike, field: str | StrField | None = None, **kwargs):
        """{}"""  # noqa: D415
        return _validate(cls=self, values=values, field=field, **kwargs)

    @doc_args(CanCurate.inspect.__doc__)
    def inspect(self, values: ListLike, field: str | StrField | None = None, **kwargs):
        """{}"""  # noqa: D415
        return _inspect(cls=self, values=values, field=field, **kwargs)

    @doc_args(CanCurate.standardize.__doc__)
    def standardize(
        self, values: Iterable, field: str | StrField | None = None, **kwargs
    ):
        """{}"""  # noqa: D415
        return _standardize(cls=self, values=values, field=field, **kwargs)


# this differs from BasicQuerySet only in .filter and .get
# QueryManager returns BasicQuerySet because it is problematic to redefine .filter and .get
# for a query set used by the default manager
class QuerySet(BasicQuerySet):
    """Sets of records returned by queries.

    Implements additional filtering capabilities.

    See Also:

        `django QuerySet <https://docs.djangoproject.com/en/4.2/ref/models/querysets/>`__

    Examples:

        >>> ULabel(name="my label").save()
        >>> queryset = ULabel.filter(name="my label")
        >>> queryset # an instance of QuerySet
    """

    def _handle_unknown_field(self, error: FieldError) -> None:
        """Suggest available fields if an unknown field was passed."""
        if "Cannot resolve keyword" in str(error):
            field = str(error).split("'")[1]
            avail_fields = self.model.__get_available_fields__()
            fields = ", ".join(sorted(avail_fields))
            raise FieldError(
                f"Unknown field '{field}'. Available fields: {fields}"
            ) from None
        raise error  # pragma: no cover

    def get(self, idlike: int | str | None = None, **expressions) -> SQLRecord:
        """Query a single record. Raises error if there are more or none."""
        is_run_input = expressions.pop("is_run_input", False)

        # artifacts_from_path and get accept only BasicQuerySet
        qs = self._to_class(BasicQuerySet, copy=True)

        if path := expressions.pop("path", None):
            from .artifact_set import ArtifactSet, artifacts_from_path

            if not isinstance(self, ArtifactSet):
                raise ValueError("Querying by path is only possible for artifacts.")
            qs = artifacts_from_path(qs, path)

        try:
            record = get(qs, idlike, **expressions)
        except ValueError as e:
            # Pass through original error for explicit id lookups
            if "Field 'id' expected a number" in str(e):
                if "id" in expressions:
                    raise
                field = next(iter(expressions))
                raise FieldError(
                    f"Invalid lookup '{expressions[field]}' for {field}. Did you mean {field}__name?"
                ) from None
            raise  # pragma: no cover
        except FieldError as e:
            self._handle_unknown_field(e)
            raise  # pragma: no cover

        if is_run_input is not False:  # might be None or True or Run
            from .artifact import Artifact, track_run_input
            from .collection import Collection

            if isinstance(record, (Artifact, Collection)):
                track_run_input(record, is_run_input)

        return record

    def filter(self, *queries, **expressions) -> QuerySet:
        """Query a set of records."""
        from lamindb.models import Artifact, Record, Run

        registry = self.model
        if not expressions.pop("_skip_filter_with_features", False) and registry in {
            Artifact,
            Run,
            Record,
        }:
            from ._feature_manager import filter_with_features

            qs = filter_with_features(self, *queries, **expressions)
        else:
            # Suggest to use __name for related fields such as id when not passed
            for field, value in expressions.items():
                if (
                    isinstance(value, str)
                    and value.strip("-").isalpha()
                    and "__" not in field
                    and hasattr(registry, field)
                ):
                    field_attr = getattr(registry, field)
                    if hasattr(field_attr, "field") and field_attr.field.related_model:
                        raise FieldError(
                            f"Invalid lookup '{value}' for {field}. Did you mean {field}__name?"
                        )
            expressions = process_expressions(self, queries, expressions)
            # need to run a query if queries or expressions are not empty
            if queries or expressions:
                try:
                    return super().filter(*queries, **expressions)
                except FieldError as e:
                    self._handle_unknown_field(e)
            qs = self
        return qs


@final
class NonInstantiableQuerySet:
    """Wrapper around QuerySet that prevents instantiation while preserving query methods."""

    def __init__(self, qs: QuerySet, registry_name: str):
        self._qs = qs
        self._name = registry_name

    def __repr__(self) -> str:
        return f"<QuerySet [{self._name}]>"

    def __call__(self, *args, **kwargs):
        raise TypeError(
            f"Cannot instantiate {self._name} from DB. "
            f"Use {self._name}.filter(), {self._name}.get(), etc. to query records."
        )

    def __getattr__(self, attr):
        return getattr(self._qs, attr)


class ModuleNamespace:
    """Namespace for accessing registries from a specific schema module.

    Args:
        query_db: Parent DB instance.
        module_name: Name of the schema module (e.g., 'bionty', 'pertdb').
    """

    def __init__(self, query_db: DB, module_name: str):
        self._query_db = query_db
        self._module_name = module_name
        self._cache: dict[str, NonInstantiableQuerySet] = {}

    def __getattr__(self, name: str) -> NonInstantiableQuerySet:
        """Access a registry class from this schema module.

        Args:
            name: Registry class name (e.g., 'Gene', 'CellType').

        Returns:
            QuerySet for the specified registry scoped to the parent instance.
        """
        if name in self._cache:
            return self._cache[name]

        try:
            schema_module = import_module(self._module_name)
            if hasattr(schema_module, name):
                model_class = getattr(schema_module, name)
                queryset = model_class.connect(self._query_db._instance)
                wrapped = NonInstantiableQuerySet(queryset, name)
                self._cache[name] = wrapped
                return wrapped
        except (ImportError, AttributeError):
            pass

        raise AttributeError(
            f"Registry '{name}' not found in lamindb. Use .bt.{name} or .pertdb.{name} for schema-specific registries."
        )

    def __dir__(self) -> list[str]:
        """Return list of available registries in this schema module."""
        base_attrs = [attr for attr in object.__dir__(self) if not attr.startswith("_")]
        try:
            schema_module = import_module(self._module_name)
            if hasattr(schema_module, "__all__"):
                registries = set()
                for class_name in schema_module.__all__:
                    model_class = getattr(schema_module, class_name, None)
                    if model_class and hasattr(model_class, "connect"):
                        registries.add(class_name)
                return sorted(set(base_attrs) | registries)
        except ImportError:
            pass
        return base_attrs


class BiontyDB(ModuleNamespace):
    """Namespace for Bionty registries (Gene, CellType, Disease, etc.)."""

    Gene: QuerySet[Gene]  # type: ignore[type-arg]
    Protein: QuerySet[Protein]  # type: ignore[type-arg]
    CellType: QuerySet[CellType]  # type: ignore[type-arg]
    Disease: QuerySet[Disease]  # type: ignore[type-arg]
    Phenotype: QuerySet[Phenotype]  # type: ignore[type-arg]
    Pathway: QuerySet[Pathway]  # type: ignore[type-arg]
    Tissue: QuerySet[Tissue]  # type: ignore[type-arg]
    CellLine: QuerySet[CellLine]  # type: ignore[type-arg]
    CellMarker: QuerySet[CellMarker]  # type: ignore[type-arg]
    Organism: QuerySet[Organism]  # type: ignore[type-arg]
    ExperimentalFactor: QuerySet[ExperimentalFactor]  # type: ignore[type-arg]
    DevelopmentalStage: QuerySet[DevelopmentalStage]  # type: ignore[type-arg]
    Ethnicity: QuerySet[Ethnicity]  # type: ignore[type-arg]


class PertdbDB(ModuleNamespace):
    """Namespace for `PertDB` registries (Biologic, Compound, etc.)."""

    Biologic: QuerySet[Biologic]  # type: ignore[type-arg]
    Compound: QuerySet[Compound]  # type: ignore[type-arg]
    CompoundPerturbation: QuerySet[CompoundPerturbation]  # type: ignore[type-arg]
    GeneticPerturbation: QuerySet[GeneticPerturbation]  # type: ignore[type-arg]
    EnvironmentalPerturbation: QuerySet[EnvironmentalPerturbation]  # type: ignore[type-arg]
    CombinationPerturbation: QuerySet[CombinationPerturbation]  # type: ignore[type-arg]
    PerturbationTarget: QuerySet[PerturbationTarget]  # type: ignore[type-arg]


class DB:
    """Query any registry of any instance.

    Args:
        instance: Instance identifier in format "account/instance".

    Examples:

        Query objects from an instance::

            db = ln.DB("laminlabs/cellxgene")

        Query artifacts and filter by `suffix`::

            db.Artifact.filter(suffix=".h5ad").to_dataframe()

        Get a single artifact by uid::

            artifact = db.Artifact.get("abcDEF123456")

        Query records and filter by name::

            db.Record.filter(name__startswith="sample").to_dataframe()

        Get a cell type object::

            t_cell = db.bionty.CellType.get(name="T cell")

        Create a lookup object to auto-complete all cell types in the database::

            cell_types = db.bionty.CellType.lookup()

        Return a `DataFrame` with additional info::

            db.Artifact.filter(
                suffix=".h5ad",
                description__contains="immune",
                size__gt=1e9,  # size > 1GB
                cell_types__name__in=["B cell", "T cell"],
            ).order_by("created_at").to_dataframe(
                include=["cell_types__name", "created_by__handle"]  # include additional info
            ).head()
    """

    Artifact: QuerySet[Artifact]  # type: ignore[type-arg]
    Collection: QuerySet[Collection]  # type: ignore[type-arg]
    Transform: QuerySet[Transform]  # type: ignore[type-arg]
    Run: QuerySet[Run]  # type: ignore[type-arg]
    User: QuerySet[User]  # type: ignore[type-arg]
    Storage: QuerySet[Storage]  # type: ignore[type-arg]
    Feature: QuerySet[Feature]  # type: ignore[type-arg]
    ULabel: QuerySet[ULabel]  # type: ignore[type-arg]
    Record: QuerySet[Record]  # type: ignore[type-arg]
    Schema: QuerySet[Schema]  # type: ignore[type-arg]
    Project: QuerySet[Project]  # type: ignore[type-arg]
    Reference: QuerySet[Reference]  # type: ignore[type-arg]
    Branch: QuerySet[Branch]  # type: ignore[type-arg]
    Space: QuerySet[Space]  # type: ignore[type-arg]

    bionty: BiontyDB
    pertdb: PertdbDB

    def __init__(self, instance: str):
        self._instance = instance
        self._cache: dict[str, NonInstantiableQuerySet | BiontyDB | PertdbDB] = {}
        self._available_registries: set[str] | None = None

        owner, instance_name = instance.split("/")
        instance_info = ln_setup._connect_instance._connect_instance(
            owner=owner, name=instance_name
        )
        self._modules = ["lamindb"] + list(instance_info.modules)

    def __getattr__(self, name: str) -> NonInstantiableQuerySet | BiontyDB | PertdbDB:
        """Access a registry class or schema namespace for this database instance.

        Args:
            name: Registry class name (e.g., 'Artifact', 'Collection') or schema namespace ('bionty', 'pertdb').

        Returns:
            QuerySet for the specified registry or schema namespace scoped to this instance.
        """
        if name in self._cache:
            return self._cache[name]

        if name == "bionty":
            if "bionty" not in self._modules:
                raise AttributeError(
                    f"Schema 'bionty' not available in instance '{self._instance}'."
                )
            if "bionty" not in self._cache:
                namespace = BiontyDB(self, "bionty")
                self._cache["bionty"] = namespace
            return self._cache["bionty"]

        if name == "pertdb":
            if "pertdb" not in self._modules:
                raise AttributeError(
                    f"Schema 'pertdb' not available in instance '{self._instance}'."
                )
            if "pertdb" not in self._cache:
                namespace = PertdbDB(self, "pertdb")  # type: ignore
                self._cache["pertdb"] = namespace
            return self._cache["pertdb"]

        try:
            lamindb_module = import_module("lamindb")
            if hasattr(lamindb_module, name):
                model_class = getattr(lamindb_module, name)
                queryset = model_class.connect(self._instance)
                wrapped = NonInstantiableQuerySet(queryset, name)
                self._cache[name] = wrapped
                return wrapped
        except (ImportError, AttributeError):
            pass

        raise AttributeError(
            f"Registry '{name}' not found in lamindb core registries. Use .bionty.{name} or .pertdb.{name} for schema-specific registries."
        )

    def __repr__(self) -> str:
        return f"DB('{self._instance}')"

    def __dir__(self) -> list[str]:
        """Return list of available registries and schema namespaces."""
        base_attrs = [attr for attr in super().__dir__() if not attr.startswith("_")]

        lamindb_registries = set()
        try:
            lamindb_module = import_module("lamindb")
            if hasattr(lamindb_module, "__all__"):
                for class_name in lamindb_module.__all__:
                    model_class = getattr(lamindb_module, class_name, None)
                    if model_class and hasattr(model_class, "connect"):
                        lamindb_registries.add(class_name)
        except ImportError:
            pass

        module_namespaces = set()
        if "bionty" in self._modules:
            module_namespaces.add("bionty")
        if "pertdb" in self._modules:
            module_namespaces.add("pertdb")

        return sorted(set(base_attrs) | lamindb_registries | module_namespaces)
