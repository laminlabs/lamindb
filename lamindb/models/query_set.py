from __future__ import annotations

import re
from collections import UserList
from collections.abc import Iterable
from collections.abc import Iterable as IterableType
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any, Generic, NamedTuple, TypeVar, Union

import pandas as pd
from django.core.exceptions import FieldError
from django.db import models
from django.db.models import F, ForeignKey, ManyToManyField, Q, Subquery
from django.db.models.fields.related import ForeignObjectRel
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args

from ..errors import DoesNotExist
from ._is_versioned import IsVersioned
from .can_curate import CanCurate, _inspect, _standardize, _validate
from .query_manager import _lookup, _search
from .sqlrecord import SQLRecord

if TYPE_CHECKING:
    from lamindb.base.types import ListLike, StrField

T = TypeVar("T")


class MultipleResultsFound(Exception):
    pass


pd.set_option("display.max_columns", 200)


# def format_and_convert_to_local_time(series: pd.Series):
#     tzinfo = datetime.now().astimezone().tzinfo
#     timedelta = tzinfo.utcoffset(datetime.now())  # type: ignore
#     return (series + timedelta).dt.strftime("%Y-%m-%d %H:%M:%S %Z")


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


def one_helper(self, does_not_exist_msg: str | None = None):
    if len(self) == 0:
        raise DoesNotExist(does_not_exist_msg)
    elif len(self) > 1:
        raise MultipleResultsFound(self)
    else:
        return self[0]


def get_backward_compat_filter_kwargs(queryset, expressions):
    from lamindb.models import (
        Artifact,
        Collection,
        Transform,
    )

    if queryset.model in {Collection, Transform}:
        name_mappings = {
            "visibility": "branch_id",
            "_branch_code": "branch_id",
        }
    elif queryset.model == Artifact:
        name_mappings = {
            "visibility": "branch_id",
            "_branch_code": "branch_id",
            "transform": "run__transform",
        }
    else:
        return expressions
    was_list = False
    if isinstance(expressions, list):
        was_list = True
        expressions = {field: True for field in expressions}
    mapped = {}
    for field, value in expressions.items():
        parts = field.split("__")
        if parts[0] in name_mappings:
            new_field = name_mappings[parts[0]] + (
                "__" + "__".join(parts[1:]) if len(parts) > 1 else ""
            )
            mapped[new_field] = value
        else:
            mapped[field] = value
    return list(mapped.keys()) if was_list else mapped


def process_expressions(queryset: QuerySet, expressions: dict) -> dict:
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

    expressions = get_backward_compat_filter_kwargs(
        queryset,
        expressions,
    )

    if issubclass(queryset.model, SQLRecord):
        # branch_id is set to 0 unless expressions contains id or uid
        if not (
            "id" in expressions
            or "uid" in expressions
            or "uid__startswith" in expressions
        ):
            branch_id = "branch_id"
            if not any(e.startswith(branch_id) for e in expressions):
                expressions[branch_id] = 1  # default branch_id
            # if branch_id is None, do not apply a filter
            # otherwise, it would mean filtering for NULL values, which doesn't make
            # sense for a non-NULLABLE column
            elif branch_id in expressions and expressions[branch_id] is None:
                expressions.pop(branch_id)
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
    registry_or_queryset: Union[type[SQLRecord], QuerySet],
    idlike: int | str | None = None,
    **expressions,
) -> SQLRecord:
    if isinstance(registry_or_queryset, QuerySet):
        qs = registry_or_queryset
        registry = qs.model
    else:
        qs = QuerySet(model=registry_or_queryset)
        registry = registry_or_queryset
    if isinstance(idlike, int):
        return super(QuerySet, qs).get(id=idlike)  # type: ignore
    elif isinstance(idlike, str):
        qs = qs.filter(uid__startswith=idlike)

        NAME_FIELD = (
            registry._name_field if hasattr(registry, "_name_field") else "name"
        )
        DOESNOTEXIST_MSG = f"No record found with uid '{idlike}'. Did you forget a keyword as in {registry.__name__}.get({NAME_FIELD}='{idlike}')?"

        if issubclass(registry, IsVersioned):
            if len(idlike) <= registry._len_stem_uid:
                return one_helper(qs.latest_version(), DOESNOTEXIST_MSG)
            else:
                return one_helper(qs, DOESNOTEXIST_MSG)
        else:
            return one_helper(qs, DOESNOTEXIST_MSG)
    else:
        assert idlike is None  # noqa: S101
        expressions = process_expressions(qs, expressions)
        # don't want branch_id here in .get(), only in .filter()
        expressions.pop("branch_id", None)
        # inject is_latest for consistency with idlike
        is_latest_was_not_in_expressions = "is_latest" not in expressions
        if issubclass(registry, IsVersioned) and is_latest_was_not_in_expressions:
            expressions["is_latest"] = True
        try:
            return registry.objects.using(qs.db).get(**expressions)
        except registry.DoesNotExist:
            # handle the case in which the is_latest injection led to a missed query
            if "is_latest" in expressions and is_latest_was_not_in_expressions:
                expressions.pop("is_latest")
                result = (
                    registry.objects.using(qs.db)
                    .filter(**expressions)
                    .order_by("-created_at")
                    .first()
                )
                if result is not None:
                    return result
            raise registry.DoesNotExist from registry.DoesNotExist


class SQLRecordList(UserList, Generic[T]):
    """Is ordered, can't be queried, but has `.df()`."""

    def __init__(self, records: Iterable[T]):
        if isinstance(records, list):
            self.data = records  # Direct assignment if already a list, no copy
        else:
            super().__init__(records)  # Let UserList handle the conversion

    def df(self) -> pd.DataFrame:
        keys = get_keys_from_df(self.data, self.data[0].__class__)
        values = [record.__dict__ for record in self.data]
        return pd.DataFrame(values, columns=keys)

    def list(
        self, field: str
    ) -> list[str]:  # meaningful to be parallel with list() in QuerySet
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
    features_input: bool | list[str],
) -> list[str]:
    exclude_field_names = ["updated_at"]
    field_names = [
        field.name
        for field in qs.model._meta.fields
        if (
            not isinstance(field, models.ForeignKey)
            and field.name not in exclude_field_names
        )
    ]
    field_names += [
        f"{field.name}_id"
        for field in qs.model._meta.fields
        if isinstance(field, models.ForeignKey)
    ]
    for field_name in [
        "version",
        "is_latest",
        "run_id",
        "created_at",
        "created_by_id",
        "updated_at",
        "_aux",
        "branch_id",
    ]:
        if field_name in field_names:
            field_names.remove(field_name)
            field_names.append(field_name)
    if field_names[0] != "uid" and "uid" in field_names:
        field_names.remove("uid")
        field_names.insert(0, "uid")
    if (
        include or features_input
    ):  # if there is features_input, reduce fields to just the first 3
        subset_field_names = field_names[:3]
        intersection = set(field_names) & set(include)
        subset_field_names += list(intersection)
        field_names = subset_field_names
    return field_names


def get_feature_annotate_kwargs(
    features: bool | list[str] | None,
) -> tuple[dict[str, Any], list[str], QuerySet]:
    from lamindb.models import (
        Artifact,
        Feature,
    )

    feature_qs = Feature.filter()
    if isinstance(features, list):
        feature_qs = feature_qs.filter(name__in=features)
        feature_names = features
    else:  # features is True -- only consider categorical features from ULabel and non-categorical features
        feature_qs = feature_qs.filter(
            Q(~Q(dtype__startswith="cat[")) | Q(dtype__startswith="cat[ULabel")
        )
        feature_names = feature_qs.list("name")
        logger.important(
            f"queried for all categorical features with dtype 'cat[ULabel...'] and non-categorical features: ({len(feature_names)}) {feature_names}"
        )
    # Get the categorical features
    cat_feature_types = {
        feature.dtype.replace("cat[", "").replace("]", "")
        for feature in feature_qs
        if feature.dtype.startswith("cat[")
    }
    # Get relationships of labels and features
    link_models_on_models = {
        getattr(
            Artifact, obj.related_name
        ).through.__get_name_with_module__(): obj.related_model.__get_name_with_module__()
        for obj in Artifact._meta.related_objects
        if obj.related_model.__get_name_with_module__() in cat_feature_types
    }
    link_models_on_models["ArtifactULabel"] = "ULabel"
    link_attributes_on_models = {
        obj.related_name: link_models_on_models[
            obj.related_model.__get_name_with_module__()
        ]
        for obj in Artifact._meta.related_objects
        if obj.related_model.__get_name_with_module__() in link_models_on_models
    }
    # Prepare Django's annotate for features
    annotate_kwargs = {}
    for link_attr, feature_type in link_attributes_on_models.items():
        annotate_kwargs[f"{link_attr}__feature__name"] = F(
            f"{link_attr}__feature__name"
        )
        field_name = (
            feature_type.split(".")[1] if "." in feature_type else feature_type
        ).lower()
        annotate_kwargs[f"{link_attr}__{field_name}__name"] = F(
            f"{link_attr}__{field_name}__name"
        )

    annotate_kwargs["_feature_values__feature__name"] = F(
        "_feature_values__feature__name"
    )
    annotate_kwargs["_feature_values__value"] = F("_feature_values__value")
    return annotate_kwargs, feature_names, feature_qs


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


def reorder_subset_columns_in_df(df: pd.DataFrame, column_order: list[str], position=3):
    valid_columns = [col for col in column_order if col in df.columns]
    all_cols = df.columns.tolist()
    remaining_cols = [col for col in all_cols if col not in valid_columns]
    new_order = remaining_cols[:position] + valid_columns + remaining_cols[position:]
    return df[new_order]


# https://lamin.ai/laminlabs/lamindata/transform/BblTiuKxsb2g0003
# https://claude.ai/chat/6ea2498c-944d-4e7a-af08-29e5ddf637d2
def reshape_annotate_result(
    df: pd.DataFrame,
    field_names: list[str],
    cols_from_include: dict[str, str] | None,
    feature_names: list[str],
    feature_qs: QuerySet | None,
) -> pd.DataFrame:
    """Reshapes tidy table to wide format.

    Args:
        field_names: List of basic fields to include in result
        df: Input dataframe with experimental data
        extra_columns: Dict specifying additional columns to process with types ('one' or 'many')
            e.g., {'ulabels__name': 'many', 'created_by__name': 'one'}
        feature_names: Feature names.
    """
    cols_from_include = cols_from_include or {}

    # initialize result with basic fields, need a copy as we're modifying it
    # will give us warnings otherwise
    result = df[field_names].copy()
    # process features if requested
    if feature_names:
        # handle feature_values
        feature_cols = ["_feature_values__feature__name", "_feature_values__value"]
        if all(col in df.columns for col in feature_cols):
            # Create two separate dataframes - one for dict values and one for non-dict values
            is_dict = df["_feature_values__value"].apply(lambda x: isinstance(x, dict))
            dict_df, non_dict_df = df[is_dict], df[~is_dict]

            # Process non-dict values using set aggregation
            non_dict_features = non_dict_df.groupby(
                ["id", "_feature_values__feature__name"]
            )["_feature_values__value"].agg(set)

            # Process dict values using first aggregation
            dict_features = dict_df.groupby(["id", "_feature_values__feature__name"])[
                "_feature_values__value"
            ].agg("first")

            # Combine the results
            combined_features = pd.concat([non_dict_features, dict_features])

            # Unstack and reset index
            feature_values = combined_features.unstack().reset_index()
            if not feature_values.empty:
                result = result.join(
                    feature_values.set_index("id"),
                    on="id",
                )

        # handle categorical features
        links_features = [
            col
            for col in df.columns
            if "feature__name" in col and col.startswith("links_")
        ]

        if links_features:
            result = process_links_features(df, result, links_features, feature_names)

        def extract_single_element(s):
            if not hasattr(s, "__len__"):  # is NaN or other scalar
                return s
            if len(s) != 1:
                # TODO: below should depend on feature._expect_many
                # logger.warning(
                #     f"expected single value because `feature._expect_many is False` but got set {len(s)} elements: {s}"
                # )
                return s
            return next(iter(s))

        for feature in feature_qs:
            if feature.name in result.columns:
                # TODO: make dependent on feature._expect_many through
                # lambda x: extract_single_element(x, feature)
                result[feature.name] = result[feature.name].apply(
                    extract_single_element
                )

        # sort columns
        result = reorder_subset_columns_in_df(result, feature_names)

    if cols_from_include:
        result = process_cols_from_include(df, result, cols_from_include)

    return result.drop_duplicates(subset=["id"])


def process_links_features(
    df: pd.DataFrame,
    result: pd.DataFrame,
    feature_cols: list[str],
    features: bool | list[str],
) -> pd.DataFrame:
    """Process links_XXX feature columns."""
    # this loops over different entities that might be linked under a feature
    for feature_col in feature_cols:
        prefix = re.match(r"links_(.+?)__feature__name", feature_col).group(1)

        value_cols = [
            col
            for col in df.columns
            if col.startswith(f"links_{prefix}__")
            and col.endswith("__name")
            and "feature__name" not in col
        ]

        if not value_cols:
            continue

        value_col = value_cols[0]
        feature_names = df[feature_col].unique()
        feature_names = feature_names[~pd.isna(feature_names)]

        # Filter features if specific ones requested
        if isinstance(features, list):
            feature_names = [f for f in feature_names if f in features]

        for feature_name in feature_names:
            mask = df[feature_col] == feature_name
            feature_values = df[mask].groupby("id")[value_col].agg(set)
            result.insert(3, feature_name, result["id"].map(feature_values))

    return result


def process_cols_from_include(
    df: pd.DataFrame, result: pd.DataFrame, extra_columns: dict[str, str]
) -> pd.DataFrame:
    """Process additional columns based on their specified types."""
    for col, col_type in extra_columns.items():
        if col not in df.columns:
            continue
        if col in result.columns:
            continue

        values = df.groupby("id")[col].agg(set if col_type == "many" else "first")
        result.insert(3, col, result["id"].map(values))

    return result


class BasicQuerySet(models.QuerySet):
    """Sets of records returned by queries.

    See Also:

        `django QuerySet <https://docs.djangoproject.com/en/stable/ref/models/querysets/>`__

    Examples:

        Any filter statement produces a query set::

            queryset = Registry.filter(name__startswith="keyword")
    """

    def __new__(cls, model=None, query=None, using=None, hints=None):
        from lamindb.models import Artifact, ArtifactSet

        # If the model is Artifact, create a new class
        # for BasicQuerySet or QuerySet that inherits from ArtifactSet.
        # This allows to add artifact specific functionality to all classes
        # inheriting from BasicQuerySet.
        # Thus all query sets of artifacts (and only of artifacts)
        # will have functions from ArtifactSet.
        if model is Artifact and not issubclass(cls, ArtifactSet):
            new_cls = type("Artifact" + cls.__name__, (cls, ArtifactSet), {})
        else:
            new_cls = cls
        return object.__new__(new_cls)

    @doc_args(SQLRecord.df.__doc__)
    def df(
        self,
        include: str | list[str] | None = None,
        features: bool | list[str] | None = None,
    ) -> pd.DataFrame:
        """{}"""  # noqa: D415
        time = datetime.now(timezone.utc)
        if include is None:
            include_input = []
        elif isinstance(include, str):
            include_input = [include]
        else:
            include_input = include
        features_input = [] if features is None else features
        include = get_backward_compat_filter_kwargs(self, include_input)
        field_names = get_basic_field_names(self, include_input, features_input)

        annotate_kwargs = {}
        feature_names: list[str] = []
        feature_qs = None
        if features:
            feature_annotate_kwargs, feature_names, feature_qs = (
                get_feature_annotate_kwargs(features)
            )
            time = logger.debug("finished feature_annotate_kwargs", time=time)
            annotate_kwargs.update(feature_annotate_kwargs)
        if include_input:
            include_input = include_input.copy()[::-1]  # type: ignore
            include_kwargs = {s: F(s) for s in include_input if s not in field_names}
            annotate_kwargs.update(include_kwargs)
        if annotate_kwargs:
            id_subquery = self.values("id")
            time = logger.debug("finished get id values", time=time)
            # for annotate, we want the queryset without filters so that joins don't affect the annotations
            query_set_without_filters = self.model.objects.filter(
                id__in=Subquery(id_subquery)
            )
            time = logger.debug("finished get query_set_without_filters", time=time)
            if self.query.order_by:
                # Apply the same ordering to the new queryset
                query_set_without_filters = query_set_without_filters.order_by(
                    *self.query.order_by
                )
                time = logger.debug("finished order by", time=time)
            queryset = query_set_without_filters.annotate(**annotate_kwargs)
            time = logger.debug("finished annotate", time=time)
        else:
            queryset = self

        df = pd.DataFrame(queryset.values(*field_names, *list(annotate_kwargs.keys())))
        if len(df) == 0:
            df = pd.DataFrame({}, columns=field_names)
            return df
        time = logger.debug("finished creating first dataframe", time=time)
        cols_from_include = analyze_lookup_cardinality(self.model, include_input)  # type: ignore
        time = logger.debug("finished analyze_lookup_cardinality", time=time)
        df_reshaped = reshape_annotate_result(
            df, field_names, cols_from_include, feature_names, feature_qs
        )
        time = logger.debug("finished reshape_annotate_result", time=time)
        pk_name = self.model._meta.pk.name
        pk_column_name = pk_name if pk_name in df.columns else f"{pk_name}_id"
        if pk_column_name in df_reshaped.columns:
            df_reshaped = df_reshaped.set_index(pk_column_name)
        time = logger.debug("finished", time=time)
        return df_reshaped

    def delete(self, *args, **kwargs):
        """Delete all records in the query set."""
        from lamindb.models import Artifact, Collection, Run, Transform

        # both Transform & Run might reference artifacts
        if self.model in {Artifact, Collection, Transform, Run}:
            for record in self:
                logger.important(f"deleting {record}")
                record.delete(*args, **kwargs)
        else:
            super().delete(*args, **kwargs)

    def list(self, field: str | None = None) -> list[SQLRecord] | list[str]:
        """Populate an (unordered) list with the results.

        Note that the order in this list is only meaningful if you ordered the underlying query set with `.order_by()`.

        Examples:
            >>> queryset.list()  # list of records
            >>> queryset.list("name")  # list of values
        """
        if field is None:
            return list(self)
        else:
            # list casting is necessary because values_list does not return a list
            return list(self.values_list(field, flat=True))

    def first(self) -> SQLRecord | None:
        """If non-empty, the first result in the query set, otherwise ``None``.

        Examples:
            >>> queryset.first()
        """
        if len(self) == 0:
            return None
        return self[0]

    def one(self) -> SQLRecord:
        """Exactly one result. Raises error if there are more or none."""
        return one_helper(self)

    def one_or_none(self) -> SQLRecord | None:
        """At most one result. Returns it if there is one, otherwise returns ``None``.

        Examples:
            >>> ULabel.filter(name="benchmark").one_or_none()
            >>> ULabel.filter(name="non existing label").one_or_none()
        """
        if len(self) == 0:
            return None
        elif len(self) == 1:
            return self[0]
        else:
            raise MultipleResultsFound(self.all())

    def latest_version(self) -> QuerySet:
        """Filter every version family by latest version."""
        if issubclass(self.model, IsVersioned):
            return self.filter(is_latest=True)
        else:
            raise ValueError("SQLRecord isn't subclass of `lamindb.core.IsVersioned`")

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
            if "_branch_code" in avail_fields:
                avail_fields.remove("_branch_code")  # backward compat
            fields = ", ".join(sorted(avail_fields))
            raise FieldError(
                f"Unknown field '{field}'. Available fields: {fields}"
            ) from None
        raise error  # pragma: no cover

    def get(self, idlike: int | str | None = None, **expressions) -> SQLRecord:
        """Query a single record. Raises error if there are more or none."""
        is_run_input = expressions.pop("is_run_input", False)

        try:
            record = get(self, idlike, **expressions)
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
            from lamindb.models.artifact import Artifact, _track_run_input
            from lamindb.models.collection import Collection

            if isinstance(record, (Artifact, Collection)):
                _track_run_input(record, is_run_input)

        return record

    def filter(self, *queries, **expressions) -> QuerySet:
        """Query a set of records."""
        # Suggest to use __name for related fields such as id when not passed
        for field, value in expressions.items():
            if (
                isinstance(value, str)
                and value.strip("-").isalpha()
                and "__" not in field
                and hasattr(self.model, field)
            ):
                field_attr = getattr(self.model, field)
                if hasattr(field_attr, "field") and field_attr.field.related_model:
                    raise FieldError(
                        f"Invalid lookup '{value}' for {field}. Did you mean {field}__name?"
                    )

        expressions = process_expressions(self, expressions)
        # need to run a query if queries or expressions are not empty
        if queries or expressions:
            try:
                return super().filter(*queries, **expressions)
            except FieldError as e:
                self._handle_unknown_field(e)
        return self
