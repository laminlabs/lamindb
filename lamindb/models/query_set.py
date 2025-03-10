from __future__ import annotations

import re
import warnings
from collections import UserList
from collections.abc import Iterable
from collections.abc import Iterable as IterableType
from typing import TYPE_CHECKING, Any, Generic, Literal, NamedTuple, TypeVar, Union

import pandas as pd
from django.core.exceptions import FieldError
from django.db import models
from django.db.models import F, ForeignKey, ManyToManyField, Subquery
from django.db.models.fields.related import ForeignObjectRel
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args

from ..core._mapped_collection import MappedCollection
from ..errors import DoesNotExist
from ._is_versioned import IsVersioned
from .can_curate import CanCurate
from .record import Record

if TYPE_CHECKING:
    import builtins

    from pyarrow.dataset import Dataset as PyArrowDataset

    from lamindb.base.types import ListLike, StrField

    from ..core.storage import UPath

T = TypeVar("T")


class MultipleResultsFound(Exception):
    pass


pd.set_option("display.max_columns", 200)


# def format_and_convert_to_local_time(series: pd.Series):
#     tzinfo = datetime.now().astimezone().tzinfo
#     timedelta = tzinfo.utcoffset(datetime.now())  # type: ignore
#     return (series + timedelta).dt.strftime("%Y-%m-%d %H:%M:%S %Z")


def get_keys_from_df(data: list, registry: Record) -> list[str]:
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


def one_helper(self):
    if len(self) == 0:
        raise DoesNotExist
    elif len(self) > 1:
        raise MultipleResultsFound(self)
    else:
        return self[0]


def get_backward_compat_filter_kwargs(queryset, expressions):
    from lamindb.models import (
        Artifact,
        Collection,
        Schema,
        Transform,
    )

    if queryset.model in {Collection, Transform}:
        name_mappings = {
            "name": "key",
            "visibility": "_branch_code",  # for convenience (and backward compat <1.0)
        }
    elif queryset.model == Artifact:
        name_mappings = {
            "n_objects": "n_files",
            "visibility": "_branch_code",  # for convenience (and backward compat <1.0)
            "transform": "run__transform",  # for convenience (and backward compat <1.0)
            "type": "kind",
            "_accessor": "otype",
        }
    elif queryset.model == Schema:
        name_mappings = {
            "registry": "itype",
        }
    else:
        return expressions
    was_list = False
    if isinstance(expressions, list):
        # make a dummy dictionary
        was_list = True
        expressions = {field: True for field in expressions}
    mapped = {}
    for field, value in expressions.items():
        parts = field.split("__")
        if parts[0] in name_mappings:
            if parts[0] not in {
                "transform",
                "visibility",
                "schemas",
                "artifacts",
            }:
                warnings.warn(
                    f"{name_mappings[parts[0]]} is deprecated, please query for {parts[0]} instead",
                    DeprecationWarning,
                    stacklevel=2,
                )
            new_field = name_mappings[parts[0]] + (
                "__" + "__".join(parts[1:]) if len(parts) > 1 else ""
            )
            mapped[new_field] = value
        else:
            mapped[field] = value
    return list(mapped.keys()) if was_list else mapped


def process_expressions(queryset: QuerySet, expressions: dict) -> dict:
    def _map_databases(value: Any, key: str, target_db: str) -> tuple[str, Any]:
        if isinstance(value, Record):
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
            if any(isinstance(v, Record) and v._state.db != target_db for v in value):
                logger.warning(
                    f"passing records from another database to query {target_db}, matching on uids"
                )
                return key.replace("__in", "__uid__in"), [
                    v.uid if isinstance(v, Record) else v for v in value
                ]
            return key, value

        return key, value

    expressions = get_backward_compat_filter_kwargs(
        queryset,
        expressions,
    )

    if issubclass(queryset.model, Record):
        # _branch_code is set to 0 unless expressions contains id or uid
        if not (
            "id" in expressions
            or "uid" in expressions
            or "uid__startswith" in expressions
        ):
            _branch_code = "_branch_code"
            if not any(e.startswith(_branch_code) for e in expressions):
                expressions[_branch_code] = 1  # default _branch_code
            # if _branch_code is None, do not apply a filter
            # otherwise, it would mean filtering for NULL values, which doesn't make
            # sense for a non-NULLABLE column
            elif _branch_code in expressions and expressions[_branch_code] is None:
                expressions.pop(_branch_code)
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
    registry_or_queryset: Union[type[Record], QuerySet],
    idlike: int | str | None = None,
    **expressions,
) -> Record:
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
        if issubclass(registry, IsVersioned):
            if len(idlike) <= registry._len_stem_uid:
                return qs.latest_version().one()
            else:
                return qs.one()
        else:
            return qs.one()
    else:
        assert idlike is None  # noqa: S101
        expressions = process_expressions(qs, expressions)
        # don't want _branch_code here in .get(), only in .filter()
        expressions.pop("_branch_code", None)
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


class RecordList(UserList, Generic[T]):
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

    def one(self) -> T:
        """Exactly one result. Throws error if there are more or none."""
        return one_helper(self)

    def save(self) -> RecordList[T]:
        """Save all records to the database."""
        from lamindb.models.save import save

        save(self)
        return self


def get_basic_field_names(
    qs: QuerySet, include: list[str], features: bool | list[str] = False
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
        "_branch_code",
    ]:
        if field_name in field_names:
            field_names.remove(field_name)
            field_names.append(field_name)
    if field_names[0] != "uid" and "uid" in field_names:
        field_names.remove("uid")
        field_names.insert(0, "uid")
    if include or features:
        subset_field_names = field_names[:4]
        intersection = set(field_names) & set(include)
        subset_field_names += list(intersection)
        field_names = subset_field_names
    return field_names


def get_feature_annotate_kwargs(show_features: bool | list[str]) -> dict[str, Any]:
    from lamindb.models import (
        Artifact,
        Feature,
    )

    features = Feature.filter()
    if isinstance(show_features, list):
        features.filter(name__in=show_features)
    # Get the categorical features
    cat_feature_types = {
        feature.dtype.replace("cat[", "").replace("]", "")
        for feature in features
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
    return annotate_kwargs


# https://claude.ai/share/16280046-6ae5-4f6a-99ac-dec01813dc3c
def analyze_lookup_cardinality(
    model_class: Record, lookup_paths: list[str] | None
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


# https://lamin.ai/laminlabs/lamindata/transform/BblTiuKxsb2g0003
# https://claude.ai/chat/6ea2498c-944d-4e7a-af08-29e5ddf637d2
def reshape_annotate_result(
    field_names: list[str],
    df: pd.DataFrame,
    extra_columns: dict[str, str] | None = None,
    features: bool | list[str] = False,
) -> pd.DataFrame:
    """Reshapes experimental data with optional feature handling.

    Parameters:
    field_names: List of basic fields to include in result
    df: Input dataframe with experimental data
    extra_columns: Dict specifying additional columns to process with types ('one' or 'many')
                  e.g., {'ulabels__name': 'many', 'created_by__name': 'one'}
    features: If False, skip feature processing. If True, process all features.
             If list of strings, only process specified features.

    Returns:
    DataFrame with reshaped data
    """
    extra_columns = extra_columns or {}

    # Initialize result with basic fields
    result = df[field_names].drop_duplicates(subset=["id"])

    # Process features if requested
    if features:
        # Handle _feature_values if columns exist
        feature_cols = ["_feature_values__feature__name", "_feature_values__value"]
        if all(col in df.columns for col in feature_cols):
            feature_values = process_feature_values(df, features)
            if not feature_values.empty:
                for col in feature_values.columns:
                    if col in result.columns:
                        continue
                    result.insert(4, col, feature_values[col])

        # Handle links features if they exist
        links_features = [
            col
            for col in df.columns
            if "feature__name" in col and col.startswith("links_")
        ]

        if links_features:
            result = process_links_features(df, result, links_features, features)

    # Process extra columns
    if extra_columns:
        result = process_extra_columns(df, result, extra_columns)

    return result


def process_feature_values(
    df: pd.DataFrame, features: bool | list[str]
) -> pd.DataFrame:
    """Process _feature_values columns."""
    feature_values = df.groupby(["id", "_feature_values__feature__name"])[
        "_feature_values__value"
    ].agg(set)

    # Filter features if specific ones requested
    if isinstance(features, list):
        feature_values = feature_values[
            feature_values.index.get_level_values(
                "_feature_values__feature__name"
            ).isin(features)
        ]

    return feature_values.unstack().reset_index()


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
            result.insert(4, feature_name, result["id"].map(feature_values))

    return result


def process_extra_columns(
    df: pd.DataFrame, result: pd.DataFrame, extra_columns: dict[str, str]
) -> pd.DataFrame:
    """Process additional columns based on their specified types."""
    for col, col_type in extra_columns.items():
        if col not in df.columns:
            continue
        if col in result.columns:
            continue

        values = df.groupby("id")[col].agg(set if col_type == "many" else "first")
        result.insert(4, col, result["id"].map(values))

    return result


class QuerySet(models.QuerySet):
    """Sets of records returned by queries.

    See Also:

        `django QuerySet <https://docs.djangoproject.com/en/4.2/ref/models/querysets/>`__

    Examples:

        >>> ULabel(name="my label").save()
        >>> queryset = ULabel.filter(name="my label")
        >>> queryset
    """

    @doc_args(Record.df.__doc__)
    def df(
        self,
        include: str | list[str] | None = None,
        features: bool | list[str] = False,
    ) -> pd.DataFrame:
        """{}"""  # noqa: D415
        if include is None:
            include = []
        elif isinstance(include, str):
            include = [include]
        include = get_backward_compat_filter_kwargs(self, include)
        field_names = get_basic_field_names(self, include, features)  # type: ignore

        annotate_kwargs = {}
        if features:
            annotate_kwargs.update(get_feature_annotate_kwargs(features))
        if include:
            include = include.copy()[::-1]  # type: ignore
            include_kwargs = {s: F(s) for s in include if s not in field_names}
            annotate_kwargs.update(include_kwargs)
        if annotate_kwargs:
            id_subquery = self.values("id")
            # for annotate, we want the queryset without filters so that joins don't affect the annotations
            query_set_without_filters = self.model.objects.filter(
                id__in=Subquery(id_subquery)
            )
            if self.query.order_by:
                # Apply the same ordering to the new queryset
                query_set_without_filters = query_set_without_filters.order_by(
                    *self.query.order_by
                )
            queryset = query_set_without_filters.annotate(**annotate_kwargs)
        else:
            queryset = self

        df = pd.DataFrame(queryset.values(*field_names, *list(annotate_kwargs.keys())))
        if len(df) == 0:
            df = pd.DataFrame({}, columns=field_names)
            return df
        extra_cols = analyze_lookup_cardinality(self.model, include)  # type: ignore
        df_reshaped = reshape_annotate_result(field_names, df, extra_cols, features)
        pk_name = self.model._meta.pk.name
        pk_column_name = pk_name if pk_name in df.columns else f"{pk_name}_id"
        if pk_column_name in df_reshaped.columns:
            df_reshaped = df_reshaped.set_index(pk_column_name)
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
            self._delete_base_class(*args, **kwargs)

    def list(self, field: str | None = None) -> list[Record]:
        """Populate a list with the results.

        Examples:
            >>> queryset.list()  # list of records
            >>> queryset.list("name")  # list of values
        """
        if field is None:
            return list(self)
        else:
            return list(self.values_list(field, flat=True))

    def first(self) -> Record | None:
        """If non-empty, the first result in the query set, otherwise ``None``.

        Examples:
            >>> queryset.first()
        """
        if len(self) == 0:
            return None
        return self[0]

    def _handle_unknown_field(self, error: FieldError) -> None:
        """Suggest available fields if an unknown field was passed."""
        if "Cannot resolve keyword" in str(error):
            field = str(error).split("'")[1]
            fields = ", ".join(
                sorted(
                    f.name
                    for f in self.model._meta.get_fields()
                    if not f.name.startswith("_")
                    and not f.name.startswith("links_")
                    and not f.name.endswith("_id")
                )
            )
            raise FieldError(
                f"Unknown field '{field}'. Available fields: {fields}"
            ) from None
        raise error  # pragma: no cover

    def get(self, idlike: int | str | None = None, **expressions) -> Record:
        """Query a single record. Raises error if there are more or none."""
        try:
            return get(self, idlike, **expressions)
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
        if len(expressions) > 0:
            try:
                return super().filter(*queries, **expressions)
            except FieldError as e:
                self._handle_unknown_field(e)
        return self

    def one(self) -> Record:
        """Exactly one result. Raises error if there are more or none."""
        return one_helper(self)

    def one_or_none(self) -> Record | None:
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
            raise ValueError("Record isn't subclass of `lamindb.core.IsVersioned`")

    def artifacts_open(self, is_run_input: bool | None = None) -> PyArrowDataset:
        from lamindb.models.artifact import Artifact, _track_run_input
        from lamindb.models.collection import _open_paths

        if self.model != Artifact:
            raise ValueError("A query set should consist of artifacts to be opened.")
        if not self.ordered:
            logger.warning(
                "this query set is unordered, consider using `.order_by()` first "
                "to avoid opening the artifacts in an arbitrary order"
            )

        artifacts: list[Artifact] = list(self)
        paths: list[UPath] = [artifact.path for artifact in artifacts]
        dataset = _open_paths(paths)
        # track only if successful
        _track_run_input(artifacts, is_run_input)
        return dataset

    # for some reason mypy considers type list[str] as QuerySet.list[str], this is why builtins.list
    def artifacts_mapped(
        self,
        layers_keys: str | builtins.list[str] | None = None,
        obs_keys: str | builtins.list[str] | None = None,
        obsm_keys: str | builtins.list[str] | None = None,
        obs_filter: dict[str, str | builtins.list[str]] | None = None,
        join: Literal["inner", "outer"] | None = "inner",
        encode_labels: bool | builtins.list[str] = True,
        unknown_label: str | dict[str, str] | None = None,
        cache_categories: bool = True,
        parallel: bool = False,
        dtype: str | None = None,
        stream: bool = False,
        is_run_input: bool | None = None,
    ) -> MappedCollection:
        from lamindb.models.artifact import Artifact, _track_run_input

        if self.model != Artifact:
            raise ValueError("A query set should consist of artifacts to be mapped.")
        if not self.ordered:
            logger.warning(
                "this query set is unordered, consider using `.order_by()` first "
                "to avoid opening the artifacts in an arbitrary order"
            )

        artifacts: list[Artifact] = []
        paths: list[UPath] = []
        for artifact in self:
            if ".h5ad" not in artifact.suffix and ".zarr" not in artifact.suffix:
                logger.warning(f"ignoring artifact with suffix {artifact.suffix}")
                continue
            elif not stream:
                paths.append(artifact.cache())
            else:
                paths.append(artifact.path)
            artifacts.append(artifact)
        ds = MappedCollection(
            paths,
            layers_keys,
            obs_keys,
            obsm_keys,
            obs_filter,
            join,
            encode_labels,
            unknown_label,
            cache_categories,
            parallel,
            dtype,
        )
        # track only if successful
        _track_run_input(artifacts, is_run_input)
        return ds


# -------------------------------------------------------------------------------------
# CanCurate
# -------------------------------------------------------------------------------------


@doc_args(Record.search.__doc__)
def search(self, string: str, **kwargs):
    """{}"""  # noqa: D415
    from .record import _search

    return _search(cls=self, string=string, **kwargs)


@doc_args(Record.lookup.__doc__)
def lookup(self, field: StrField | None = None, **kwargs) -> NamedTuple:
    """{}"""  # noqa: D415
    from .record import _lookup

    return _lookup(cls=self, field=field, **kwargs)


@doc_args(CanCurate.validate.__doc__)
def validate(self, values: ListLike, field: str | StrField | None = None, **kwargs):
    """{}"""  # noqa: D415
    from .can_curate import _validate

    return _validate(cls=self, values=values, field=field, **kwargs)


@doc_args(CanCurate.inspect.__doc__)
def inspect(self, values: ListLike, field: str | StrField | None = None, **kwargs):
    """{}"""  # noqa: D415
    from .can_curate import _inspect

    return _inspect(cls=self, values=values, field=field, **kwargs)


@doc_args(CanCurate.standardize.__doc__)
def standardize(self, values: Iterable, field: str | StrField | None = None, **kwargs):
    """{}"""  # noqa: D415
    from .can_curate import _standardize

    return _standardize(cls=self, values=values, field=field, **kwargs)


models.QuerySet.df = QuerySet.df
models.QuerySet.list = QuerySet.list
models.QuerySet.first = QuerySet.first
models.QuerySet.one = QuerySet.one
models.QuerySet.one_or_none = QuerySet.one_or_none
models.QuerySet.latest_version = QuerySet.latest_version
models.QuerySet.search = search
models.QuerySet.lookup = lookup
models.QuerySet.validate = validate
models.QuerySet.inspect = inspect
models.QuerySet.standardize = standardize
models.QuerySet._delete_base_class = models.QuerySet.delete
models.QuerySet.delete = QuerySet.delete
