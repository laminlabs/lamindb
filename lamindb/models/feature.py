from __future__ import annotations

import importlib
from typing import TYPE_CHECKING, Any, get_args, overload

import numpy as np
import pandas as pd
from django.db import models
from django.db.models import CASCADE, PROTECT
from django.db.models.query_utils import DeferredAttribute
from django.db.utils import IntegrityError
from lamin_utils import logger
from lamindb_setup._init_instance import get_schema_module_name
from lamindb_setup.core.hashing import HASH_LENGTH, hash_dict, hash_string
from lamindb_setup.errors import (
    MODULE_WASNT_CONFIGURED_MESSAGE_TEMPLATE,
    ModuleWasntConfigured,
)
from pandas.api.types import CategoricalDtype, is_string_dtype
from pandas.core.dtypes.base import ExtensionDtype

from lamindb.base.fields import (
    BooleanField,
    CharField,
    ForeignKey,
    JSONField,
    TextField,
)
from lamindb.base.types import Dtype, FieldAttr
from lamindb.errors import DoesNotExist, FieldValidationError, ValidationError

from ..base.ids import base62_12
from ._relations import dict_module_name_to_model_name
from .can_curate import CanCurate
from .query_set import SQLRecordList
from .run import (
    TracksRun,
    TracksUpdates,
)
from .sqlrecord import BaseSQLRecord, Registry, SQLRecord, _get_record_kwargs

if TYPE_CHECKING:
    from collections.abc import Iterable

    from .schema import Schema

FEATURE_DTYPES = set(get_args(Dtype))


def parse_dtype(dtype_str: str, is_param: bool = False) -> list[dict[str, Any]]:
    """Parses feature data type string into a structured list of components."""
    from .artifact import Artifact

    allowed_dtypes = FEATURE_DTYPES
    if is_param:
        allowed_dtypes.add("dict")

    # Handle list[...] types
    if dtype_str.startswith("list[") and dtype_str.endswith("]"):
        inner_dtype_str = dtype_str[5:-1]  # Remove "list[" and "]"
        # Recursively parse the inner type
        inner_result = parse_dtype(inner_dtype_str, is_param)
        # Add "list": True to each component
        for component in inner_result:
            if isinstance(component, dict):
                component["list"] = True  # type: ignore
        return inner_result

    is_composed_cat = dtype_str.startswith("cat[") and dtype_str.endswith("]")
    result = []
    if is_composed_cat:
        related_registries = dict_module_name_to_model_name(Artifact)
        registries_str = dtype_str.replace("cat[", "")[:-1]  # strip last ]
        if registries_str != "":
            registry_str_list = registries_str.split("|")
            for cat_single_dtype_str in registry_str_list:
                single_result = parse_cat_dtype(
                    cat_single_dtype_str, related_registries
                )
                result.append(single_result)
    elif dtype_str not in allowed_dtypes:
        raise ValueError(
            f"dtype is '{dtype_str}' but has to be one of {FEATURE_DTYPES}!"
        )
    return result


def parse_cat_dtype(
    dtype_str: str,
    related_registries: dict[str, SQLRecord] | None = None,
    is_itype: bool = False,
) -> dict[str, Any]:
    """Parses a categorical dtype string into its components (registry, field, subtypes)."""
    from .artifact import Artifact

    assert isinstance(dtype_str, str)  # noqa: S101
    if related_registries is None:
        related_registries = dict_module_name_to_model_name(Artifact)

    # Parse the string considering nested brackets
    parsed = parse_nested_brackets(dtype_str)

    registry_str = parsed["registry"]
    sub_type_str = parsed["subtype"]
    field_str = parsed["field"]

    if not is_itype:
        if registry_str not in related_registries:
            raise ValidationError(
                f"'{registry_str}' is an invalid dtype, has to be registry, e.g. ULabel or bionty.CellType"
            )
        registry = related_registries[registry_str]
    else:
        if "." in registry_str:
            registry_str_split = registry_str.split(".")
            assert len(registry_str_split) == 2, registry_str  # noqa: S101
            module_name_attempt, class_name = registry_str_split
            module_name = get_schema_module_name(
                module_name_attempt, raise_import_error=False
            )
            if module_name is None:
                raise ModuleWasntConfigured(
                    MODULE_WASNT_CONFIGURED_MESSAGE_TEMPLATE.format(module_name_attempt)
                )
        else:
            module_name, class_name = "lamindb", registry_str
        module = importlib.import_module(module_name)
        registry = getattr(module, class_name)

    if sub_type_str != "":
        pass
        # validate that the subtype is a record in the registry with is_type = True
    if field_str != "":
        pass
        # validate that field_str is an actual field of the module
    else:
        field_str = registry._name_field if hasattr(registry, "_name_field") else "name"

    result = {
        "registry": registry,  # should be typed as CanCurate
        "registry_str": registry_str,
        "subtype_str": sub_type_str,
        "field_str": field_str,
        "field": getattr(registry, field_str),
    }

    # Add nested subtype information if present
    if parsed.get("nested_subtypes"):
        result["nested_subtypes"] = parsed["nested_subtypes"]

    return result


def parse_nested_brackets(dtype_str: str) -> dict[str, str]:
    """Parse dtype string with potentially nested brackets.

    Examples:
        "A" -> {"registry": "A", "subtype": "", "field": ""}
        "A.field" -> {"registry": "A", "subtype": "", "field": "field"}
        "A[B]" -> {"registry": "A", "subtype": "B", "field": ""}
        "A[B].field" -> {"registry": "A", "subtype": "B", "field": "field"}
        "A[B[C]]" -> {"registry": "A", "subtype": "B[C]", "field": "", "nested_subtypes": ["B", "C"]}
        "A[B[C]].field" -> {"registry": "A", "subtype": "B[C]", "field": "field", "nested_subtypes": ["B", "C"]}

    Args:
        dtype_str: The dtype string to parse

    Returns:
        Dictionary with parsed components
    """
    if "[" not in dtype_str:
        # No brackets - handle simple cases like "A" or "A.field"
        if "." in dtype_str:
            parts = dtype_str.split(".")
            if len(parts) == 2 and parts[1][0].isupper():
                # bionty.CellType
                return {"registry": dtype_str, "subtype": "", "field": ""}
            elif len(parts) == 3:
                # bionty.CellType.name
                return {
                    "registry": f"{parts[0]}.{parts[1]}",
                    "subtype": "",
                    "field": parts[2],
                }
            else:
                # ULabel.name
                return {"registry": parts[0], "subtype": "", "field": parts[1]}
        else:
            # Simple registry name
            return {"registry": dtype_str, "subtype": "", "field": ""}

    # Find the first opening bracket
    first_bracket = dtype_str.index("[")
    # Handle case where registry_part contains a field (e.g., "bionty.Gene.ensembl_gene_id[filters]")
    registry_and_field = dtype_str[:first_bracket]
    if "." in registry_and_field:
        parts = registry_and_field.split(".")
        if len(parts) == 3:
            registry_part = f"{parts[0]}.{parts[1]}"
            pre_bracket_field = parts[2]
        else:
            registry_part = registry_and_field
            pre_bracket_field = ""
    else:
        registry_part = registry_and_field
        pre_bracket_field = ""

    # Find the matching closing bracket for the first opening bracket
    bracket_count = 0
    closing_bracket_pos = -1

    for i in range(first_bracket, len(dtype_str)):
        if dtype_str[i] == "[":
            bracket_count += 1
        elif dtype_str[i] == "]":
            bracket_count -= 1
            if bracket_count == 0:
                closing_bracket_pos = i
                break

    if closing_bracket_pos == -1:
        raise ValueError(f"Unmatched brackets in dtype string: {dtype_str}")

    # Extract subtype (everything between first [ and matching ])
    subtype_part = dtype_str[first_bracket + 1 : closing_bracket_pos]

    # Check for field after the closing bracket
    field_part = ""
    remainder = dtype_str[closing_bracket_pos + 1 :]
    if remainder.startswith("."):
        field_part = remainder[1:]  # Remove the dot

    # Use pre_bracket_field if no post_bracket field
    if not field_part and pre_bracket_field:
        field_part = pre_bracket_field

    result = {"registry": registry_part, "subtype": subtype_part, "field": field_part}

    # If subtype contains brackets, extract nested subtypes for reference
    if "[" in subtype_part:
        nested_subtypes = extract_nested_subtypes(subtype_part)
        if nested_subtypes:
            result["nested_subtypes"] = nested_subtypes  # type: ignore

    return result


def extract_nested_subtypes(subtype_str: str) -> list[str]:
    """Extract all nested subtype levels from a nested subtype string.

    Examples:
        "B[C]" -> ["B", "C"]
        "B[C[D]]" -> ["B", "C", "D"]
        "B[C[D[E]]]" -> ["B", "C", "D", "E"]

    Args:
        subtype_str: The subtype string with potential nesting

    Returns:
        List of subtype levels from outermost to innermost
    """
    subtypes = []
    current = subtype_str

    while "[" in current:
        # Find the first part before the bracket
        bracket_pos = current.index("[")
        subtypes.append(current[:bracket_pos])

        # Find the matching closing bracket
        bracket_count = 0
        closing_pos = -1

        for i in range(bracket_pos, len(current)):
            if current[i] == "[":
                bracket_count += 1
            elif current[i] == "]":
                bracket_count -= 1
                if bracket_count == 0:
                    closing_pos = i
                    break

        if closing_pos == -1:
            break

        # Move to the content inside the brackets
        current = current[bracket_pos + 1 : closing_pos]

    # Add the final innermost subtype
    if current:
        subtypes.append(current)

    return subtypes


def serialize_dtype(
    dtype: Registry
    | SQLRecord
    | FieldAttr
    | list[SQLRecord]
    | list[Registry]
    | list[str]
    | list[float]
    | str
    | type,
    is_itype: bool = False,
) -> str:
    """Converts a data type object into its string representation."""
    from .record import Record
    from .ulabel import ULabel

    # Handle generic types like list[str], list[Registry], etc.
    if hasattr(dtype, "__origin__") and dtype.__origin__ is list:
        # Get the inner type from list[T]
        inner_type = dtype.__args__[0] if dtype.__args__ else None  # type: ignore
        if inner_type is not None:
            # Recursively serialize the inner type
            inner_dtype_str = serialize_dtype(inner_type, is_itype=is_itype)
            return f"list[{inner_dtype_str}]"

    if (
        not isinstance(dtype, list)
        and hasattr(dtype, "__name__")
        and dtype.__name__ in FEATURE_DTYPES
    ):
        dtype_str = dtype.__name__
    elif dtype is dict:
        dtype_str = "dict"
    elif is_itype and isinstance(dtype, str):
        if dtype not in "Feature":
            parse_cat_dtype(
                dtype_str=dtype, is_itype=True
            )  # throws an error if invalid
        dtype_str = dtype
    elif isinstance(dtype, (ExtensionDtype, np.dtype)):
        dtype_str = serialize_pandas_dtype(dtype)
    else:
        error_message = "dtype has to be a registry, a ulabel subtype, a registry field, or a list of registries or fields, not {}"
        if isinstance(dtype, (Registry, DeferredAttribute, ULabel, Record)):
            dtype = [dtype]
        elif not isinstance(dtype, list):
            raise ValueError(error_message.format(dtype))
        dtype_str = ""
        for one_dtype in dtype:
            if not isinstance(one_dtype, (Registry, DeferredAttribute, ULabel, Record)):
                raise ValueError(error_message.format(one_dtype))
            if isinstance(one_dtype, Registry):
                dtype_str += one_dtype.__get_name_with_module__() + "|"
            elif isinstance(one_dtype, (ULabel, Record)):
                assert one_dtype.is_type, (  # noqa: S101
                    f"ulabel has to be a type if acting as dtype, {one_dtype} has `is_type` False"
                )
                if isinstance(one_dtype, ULabel):
                    dtype_str += f"ULabel[{one_dtype.name}]"
                else:
                    dtype_str += f"Record[{one_dtype.name}]"
            else:
                name = one_dtype.field.name
                field_ext = f".{name}" if name != "name" else ""
                dtype_str += (
                    one_dtype.field.model.__get_name_with_module__() + field_ext + "|"
                )
        dtype_str = dtype_str.rstrip("|")
        if not is_itype:
            dtype_str = f"cat[{dtype_str}]"
    return dtype_str


def serialize_pandas_dtype(pandas_dtype: ExtensionDtype) -> str:
    """Convert pandas ExtensionDtype to simplified string representation."""
    if is_string_dtype(pandas_dtype):
        if not isinstance(pandas_dtype, CategoricalDtype):
            dtype = "str"
        else:
            dtype = "cat"
    # there are string-like categoricals and "pure" categoricals (pd.Categorical)
    elif isinstance(pandas_dtype, CategoricalDtype):
        dtype = "cat"
    else:
        # strip precision qualifiers
        dtype = "".join(dt for dt in pandas_dtype.name if not dt.isdigit())
        if dtype == "uint":
            dtype = "int"
    if dtype.startswith("datetime"):
        dtype = dtype.split("[")[0]
    assert dtype in FEATURE_DTYPES  # noqa: S101
    return dtype


def parse_filter_string(filter_str: str) -> dict[str, tuple[str, str | None, str]]:
    """Parse comma-separated Django filter expressions into structured components.

    Args:
        filter_str: Comma-separated filters like 'name=value, relation__field=value'

    Returns:
        Dict mapping original filter key to (relation_name, field_name, value) tuple.
        For direct fields: field_name is None.
        For relations: field_name contains the lookup field.
    """
    filters = {}

    filter_parts = [part.strip() for part in filter_str.split(",")]
    for part in filter_parts:
        if "=" not in part:
            raise ValueError(f"Invalid filter expression: '{part}' (missing '=' sign)")

        key, value = part.split("=", 1)
        key = key.strip()
        value = value.strip().strip("'\"")

        if not key:
            raise ValueError(f"Invalid filter expression: '{part}' (empty key)")
        if not value:
            raise ValueError(f"Invalid filter expression: '{part}' (empty value)")

        if "__" in key:
            relation_name, field_name = key.split("__", 1)
            filters[key] = (relation_name, field_name, value)
        else:
            filters[key] = (key, None, value)

    return filters


def resolve_relation_filters(
    parsed_filters: dict[str, tuple[str, str | None, str]], registry: SQLRecord
) -> dict[str, str | SQLRecord]:
    """Resolve relation filters actual model objects.

    Args:
        parsed_filters: Django filters like output from :func:`lamindb.models.feature.parse_filter_string`
        registry: Model class to resolve relationships against

    Returns:
        Dict with resolved objects for successful relations, original values for direct fields and failed resolutions.
    """
    resolved = {}

    for filter_key, (relation_name, field_name, value) in parsed_filters.items():
        if field_name is not None:  # relation filter
            if hasattr(registry, relation_name):
                relation_field = getattr(registry, relation_name)
                if (
                    hasattr(relation_field, "field")
                    and relation_field.field.is_relation
                ):
                    try:
                        related_model = relation_field.field.related_model
                        related_obj = related_model.get(**{field_name: value})
                        resolved[relation_name] = related_obj
                        continue
                    except (DoesNotExist, AttributeError):
                        pass  # Fall back to original filter
        resolved[filter_key] = value

    return resolved


def process_init_feature_param(args, kwargs, is_param: bool = False):
    # now we proceed with the user-facing constructor
    if len(args) != 0:
        raise ValueError("Only keyword args allowed")
    name: str = kwargs.pop("name", None)
    dtype: type | str | None = kwargs.pop("dtype", None)
    is_type: bool = kwargs.pop("is_type", None)
    type_: Feature | str | None = kwargs.pop("type", None)
    description: str | None = kwargs.pop("description", None)
    branch = kwargs.pop("branch", None)
    branch_id = kwargs.pop("branch_id", 1)
    space = kwargs.pop("space", None)
    space_id = kwargs.pop("space_id", 1)
    _skip_validation = kwargs.pop("_skip_validation", False)
    if kwargs:
        valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Feature)])
        raise FieldValidationError(f"Only {valid_keywords} are valid keyword arguments")
    kwargs["name"] = name
    kwargs["type"] = type_
    kwargs["is_type"] = is_type
    kwargs["branch"] = branch
    kwargs["branch_id"] = branch_id
    kwargs["space"] = space
    kwargs["space_id"] = space_id
    kwargs["_skip_validation"] = _skip_validation
    if not is_param:
        kwargs["description"] = description
    # cast dtype
    if dtype is None and not is_type:
        raise ValidationError(
            f"Please pass dtype, one of {FEATURE_DTYPES} or a composed categorical dtype"
        )
    dtype_str = None
    if dtype is not None:
        if not isinstance(dtype, str):
            dtype_str = serialize_dtype(dtype)
        else:
            dtype_str = dtype
            parse_dtype(dtype_str, is_param=is_param)
        kwargs["dtype"] = dtype_str
    return kwargs


class Feature(SQLRecord, CanCurate, TracksRun, TracksUpdates):
    """Variables, such as dataframe columns or run parameters.

    A feature often represents a dimension of a dataset, such as a column in a `DataFrame`.
    The `Feature` registry organizes metadata of features.

    The `Feature` registry helps you organize and query datasets based on their
    features and corresponding label annotations. For instance, when working
    with a "T cell" label, it could be measured through different features
    such as `"cell_type_by_expert"` where an expert manually classified the
    cell, or `"cell_type_by_model"` where a computational model made the classification.

    The two most important metadata of a feature are its `name` and the `dtype`.
    In addition to typical data types, LaminDB has a `"num"` `dtype` to
    concisely denote the union of all numerical types.

    Args:
        name: `str` Name of the feature, typically a column name.
        dtype: `Dtype | Registry | list[Registry] | FieldAttr` See :class:`~lamindb.base.types.Dtype`.
            For categorical types, you can define to which registry values are
            restricted, e.g., `ULabel` or `[ULabel, bionty.CellType]`.
        unit: `str | None = None` Unit of measure, ideally SI (`"m"`, `"s"`, `"kg"`, etc.) or `"normalized"` etc.
        description: `str | None = None` A description.
        synonyms: `str | None = None` Bar-separated synonyms.
        nullable: `bool = True` Whether the feature can have null-like values (`None`, `pd.NA`, `NaN`, etc.), see :attr:`~lamindb.Feature.nullable`.
        default_value: `Any | None = None` Default value for the feature.
        coerce_dtype: `bool = False` When True, attempts to coerce values to the specified dtype
            during validation, see :attr:`~lamindb.Feature.coerce_dtype`.
        cat_filters: `dict[str, str] | None = None` Subset a registry by additional filters to define valid categories.

    Note:

        For more control, you can use :mod:`bionty` registries to manage simple
        biological entities like genes, proteins & cell markers. Or you define
        custom registries to manage high-level derived features like gene sets.

    See Also:
        :meth:`~lamindb.Feature.from_df`
            Create feature records from DataFrame.
        :attr:`~lamindb.Artifact.features`
            Feature manager of an artifact or collection.
        :class:`~lamindb.ULabel`
            Universal labels.
        :class:`~lamindb.Schema`
            Feature sets.

    Example:

        A simple `"str"` feature.::

            ln.Feature(name="sample_note", dtype=str).save()

        A dtype `"cat[ULabel]"` can be more easily passed as below.::

            ln.Feature(name="project", dtype=ln.ULabel).save()

        A dtype `"cat[ULabel|bionty.CellType]"` can be more easily passed as below.::

            ln.Feature(
                name="cell_type",
                dtype=[ln.ULabel, bt.CellType],
            ).save()

        A multivalue feature with a list of cell types.::

            ln.Feature(
                name="cell_types",
                dtype=list[bt.CellType],  # or list[str] for a list of strings
            ).save()

        A path feature.::

            ln.Feature(
                name="image_path",
                dtype="path",   # will be validated as `str`
            ).save()

    Hint:

        *Features* and *labels* denote two ways of using entities to organize data:

        1. A feature qualifies *what* is measured, i.e., a numerical or categorical random variable
        2. A label *is* a measured value, i.e., a category

        Consider annotating a dataset by that it measured expression of 30k
        genes: genes relate to the dataset as feature identifiers through a
        feature set with 30k members. Now consider annotating the artifact by
        whether that it measured the knock-out of 3 genes: here, the 3 genes act
        as labels of the dataset.

        Re-shaping data can introduce ambiguity among features & labels. If this
        happened, ask yourself what the joint measurement was: a feature
        qualifies variables in a joint measurement. The canonical data matrix
        lists jointly measured variables in the columns.
    """

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "name"
    _aux_fields: dict[str, tuple[str, type]] = {
        "0": ("default_value", Any),  # type: ignore
        "1": ("nullable", bool),
        "2": ("coerce_dtype", bool),
    }

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=12, default=base62_12
    )
    """Universal id, valid across DB instances."""
    name: str = CharField(max_length=150, db_index=True)
    """Name of feature."""
    dtype: Dtype | None = CharField(db_index=True, null=True)
    """Data type (:class:`~lamindb.base.types.Dtype`)."""
    type: Feature | None = ForeignKey(
        "self", PROTECT, null=True, related_name="features"
    )
    """Type of feature (e.g., 'Readout', 'Metric', 'Metadata', 'ExpertAnnotation', 'ModelPrediction').

    Allows to group features by type, e.g., all read outs, all metrics, etc.
    """
    features: Feature
    """Features of this type (can only be non-empty if `is_type` is `True`)."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    unit: str | None = CharField(max_length=30, db_index=True, null=True)
    """Unit of measure, ideally SI (`m`, `s`, `kg`, etc.) or 'normalized' etc. (optional)."""
    description: str | None = CharField(db_index=True, null=True)
    """A description."""
    array_rank: int = models.SmallIntegerField(default=0, db_index=True)
    """Rank of feature.

    Number of indices of the array: 0 for scalar, 1 for vector, 2 for matrix.

    Is called `.ndim` in `numpy` and `pytorch` but shouldn't be confused with
    the dimension of the feature space.
    """
    array_size: int = models.IntegerField(default=0, db_index=True)
    """Number of elements of the feature.

    Total number of elements (product of shape components) of the array.

    - A number or string (a scalar): 1
    - A 50-dimensional embedding: 50
    - A 25 x 25 image: 625
    """
    array_shape: list[int] | None = JSONField(default=None, db_default=None, null=True)
    """Shape of the feature.

    - A number or string (a scalar): [1]
    - A 50-dimensional embedding: [50]
    - A 25 x 25 image: [25, 25]

    Is stored as a list rather than a tuple because it's serialized as JSON.
    """
    proxy_dtype: Dtype | None = CharField(default=None, null=True)
    """Proxy data type.

    If the feature is an image it's often stored via a path to the image file. Hence, while the dtype might be
    image with a certain shape, the proxy dtype would be str.
    """
    synonyms: str | None = TextField(null=True)
    """Bar-separated (|) synonyms (optional)."""
    # we define the below ManyToMany on the feature model because it parallels
    # how other registries (like Gene, Protein, etc.) relate to Schema
    # it makes the API more consistent
    schemas: Schema = models.ManyToManyField(
        "Schema", through="SchemaFeature", related_name="features"
    )
    """Feature sets linked to this feature."""
    _expect_many: bool = models.BooleanField(default=None, db_default=None, null=True)
    """Indicates whether values for this feature are expected to occur a single or multiple times for an artifact (default `None`).

    - if it's `True` (default), the values come from an observation-level aggregation and a dtype of `datetime` on the observation-level means `set[datetime]` on the artifact-level
    - if it's `False` it's an artifact-level value and datetime means datetime; this is an edge case because an arbitrary artifact would always be a set of arbitrary measurements that would need to be aggregated ("one just happens to measure a single cell line in that artifact")
    """
    _curation: dict[str, Any] = JSONField(default=None, db_default=None, null=True)
    # backward fields
    values: FeatureValue
    """Values for this feature."""

    @overload
    def __init__(
        self,
        name: str,
        dtype: Dtype | Registry | list[Registry] | FieldAttr,
        type: Feature | None = None,
        is_type: bool = False,
        unit: str | None = None,
        description: str | None = None,
        synonyms: str | None = None,
        nullable: bool = True,
        default_value: str | None = None,
        coerce_dtype: bool = False,
        cat_filters: dict[str, str] | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        default_value = kwargs.pop("default_value", None)
        nullable = kwargs.pop("nullable", True)  # default value of nullable
        cat_filters = kwargs.pop("cat_filters", None)
        coerce_dtype = kwargs.pop("coerce_dtype", False)
        kwargs = process_init_feature_param(args, kwargs)
        super().__init__(*args, **kwargs)
        self.default_value = default_value
        self.nullable = nullable
        self.coerce_dtype = coerce_dtype
        dtype_str = kwargs.pop("dtype", None)
        if cat_filters:
            if "|" in dtype_str:
                raise ValidationError(
                    f"cat_filters are incompatible with union dtypes: '{dtype_str}'"
                )
            if "]]" in dtype_str:
                raise ValidationError(
                    f"cat_filters are incompatible with nested dtypes: '{dtype_str}'"
                )

            # Validate filter values and SQLRecord attributes
            for filter_key, filter_value in cat_filters.items():
                if not filter_value or (
                    isinstance(filter_value, str) and not filter_value.strip()
                ):
                    raise ValidationError(f"Empty value in filter {filter_key}")
                # Check SQLRecord attributes for relation lookups
                if isinstance(filter_value, SQLRecord) and "__" in filter_key:
                    field_name = filter_key.split("__", 1)[1]
                    if not hasattr(filter_value, field_name):
                        raise ValidationError(
                            f"SQLRecord {filter_value.__class__.__name__} has no attribute '{field_name}' in filter {filter_key}"
                        )

            # If a SQLRecord is passed, we access its uid to apply a standard filter
            cat_filters = {
                f"{key}__uid"
                if (
                    is_sqlrecord := isinstance(filter, SQLRecord)
                    and hasattr(filter, "uid")
                )
                else key: filter.uid if is_sqlrecord else filter
                for key, filter in cat_filters.items()
            }

            fill_in = ", ".join(
                f"{key}='{value}'" for (key, value) in cat_filters.items()
            )
            dtype_str = dtype_str.replace("]", f"[{fill_in}]]")
            self.dtype = dtype_str
        if not self._state.adding:
            if not (
                self.dtype.startswith("cat")
                if dtype_str == "cat"
                else dtype_str.startswith("cat")
                if self.dtype == "cat"
                else self.dtype == dtype_str
            ):
                raise ValidationError(
                    f"Feature {self.name} already exists with dtype {self.dtype}, you passed {dtype_str}"
                )

    @classmethod
    def from_df(cls, df: pd.DataFrame, field: FieldAttr | None = None) -> SQLRecordList:
        """Create Feature records for columns."""
        field = Feature.name if field is None else field
        registry = field.field.model  # type: ignore
        if registry != Feature:
            raise ValueError("field must be a Feature FieldAttr!")
        categoricals = categoricals_from_df(df)
        dtypes = {}
        for name, col in df.items():
            if name in categoricals:
                dtypes[name] = "cat"
            else:
                dtypes[name] = serialize_pandas_dtype(col.dtype)
        with logger.mute():  # silence the warning "loaded record with exact same name "
            features = [
                Feature(name=name, dtype=dtype) for name, dtype in dtypes.items()
            ]  # type: ignore
        assert len(features) == len(df.columns)  # noqa: S101
        return SQLRecordList(features)

    def save(self, *args, **kwargs) -> Feature:
        """Save."""
        super().save(*args, **kwargs)
        return self

    def with_config(self, optional: bool | None = None) -> tuple[Feature, dict]:
        """Pass addtional configurations to the schema."""
        if optional is not None:
            return self, {"optional": optional}
        return self, {}

    @property
    def default_value(self) -> Any:
        """A default value that overwrites missing values (default `None`).

        This takes effect when you call `Curator.standardize()`.

        If `default_value = None`, missing values like `pd.NA` or `np.nan` are kept.
        """
        if self._aux is not None and "af" in self._aux and "0" in self._aux["af"]:  # type: ignore
            return self._aux["af"]["0"]  # type: ignore
        else:
            return None

    @default_value.setter
    def default_value(self, value: str | None) -> None:
        self._aux = self._aux or {}
        self._aux.setdefault("af", {})["0"] = value

    @property
    def nullable(self) -> bool:
        """Indicates whether the feature can have nullable values (default `True`).

        Example::

            import lamindb as ln
            import pandas as pd

            disease = ln.Feature(name="disease", dtype=ln.ULabel, nullable=False).save()
            schema = ln.Schema(features=[disease]).save()
            dataset = {"disease": pd.Categorical([pd.NA, "asthma"])}
            df = pd.DataFrame(dataset)
            curator = ln.curators.DataFrameCurator(df, schema)
            try:
                curator.validate()
            except ln.errors.ValidationError as e:
                assert str(e).startswith("non-nullable series 'disease' contains null values")

        """
        if self._aux is not None and "af" in self._aux and "1" in self._aux["af"]:
            value = self._aux["af"]["1"]
            return True if value is None else value
        else:
            return True

    @nullable.setter
    def nullable(self, value: bool) -> None:
        assert isinstance(value, bool), value  # noqa: S101
        self._aux = self._aux or {}
        self._aux.setdefault("af", {})["1"] = value

    @property
    def coerce_dtype(self) -> bool:
        """Whether dtypes should be coerced during validation.

        For example, a `objects`-dtyped pandas column can be coerced to `categorical` and would pass validation if this is true.
        """
        if self._aux is not None and "af" in self._aux and "2" in self._aux["af"]:  # type: ignore
            return self._aux["af"]["2"]  # type: ignore
        else:
            return False

    @coerce_dtype.setter
    def coerce_dtype(self, value: bool) -> None:
        self._aux = self._aux or {}
        self._aux.setdefault("af", {})["2"] = value

    # we'll enable this later
    # @property
    # def observational_unit(self) -> Literal["Artifact", "Observation"]:
    #     """Default observational unit on which the feature is measured.

    #     Currently, we only make a distinction between artifact-level and observation-level features.

    #     For example, a feature `"ml_split"` that stores `"test"` & `"train"` labels is typically defined on the artifact level.
    #     When accessing `artifact.features.get_values(["ml_split"])`, you expect a single value, either `"test"` or `"train"`.

    #     However, when accessing an artifact annotation with a feature that's defined on the observation-level, say `"cell_type"`, you expect a set of values. So,
    #     `artifact.features.get_values(["cell_type_from_expert"])` should return a set: `{"T cell", "B cell"}`.

    #     The value of `observational_unit` is currently auto-managed: if using `artifact.featueres.add_values()`,
    #     it will be set to `Artifact`. In a curator, the value depends on whether it's an artifact- or observation-level slot
    #     (e.g. `.uns` is artifact-level in `AnnData` whereas `.obs` is observation-level).

    #     Note: This attribute might in the future be used to distinguish different types of observational units (e.g. single cells vs. physical samples vs. study subjects etc.).
    #     """
    #     if self._expect_many:
    #         return "Observation"  # this here might be replaced with the specific observational unit
    #     else:
    #         return "Artifact"


class FeatureValue(SQLRecord, TracksRun):
    """Non-categorical features values.

    Categorical feature values are stored in their respective registries:
    :class:`~lamindb.ULabel`, :class:`~bionty.CellType`, etc.

    Unlike for ULabel, in `FeatureValue`, values are grouped by features and
    not by an ontological hierarchy.
    """

    # we do not have a unique constraint on feature & value because it leads to hashing errors
    # for large dictionaries: https://lamin.ai/laminlabs/lamindata/transform/jgTrkoeuxAfs0000
    # we do not hash values because we have `get_or_create` logic all over the place
    # and also for checking whether the (feature, value) combination exists
    # there does not seem an issue with querying for a dict-like value
    # https://lamin.ai/laminlabs/lamindata/transform/jgTrkoeuxAfs0001

    _name_field: str = "value"

    feature: Feature | None = ForeignKey(
        Feature, CASCADE, null=True, related_name="values", default=None
    )
    """The dimension metadata."""
    value: Any = models.JSONField()
    """The JSON-like value."""
    hash: str = CharField(max_length=HASH_LENGTH, null=True, db_index=True)
    """Value hash."""

    class Meta(BaseSQLRecord.Meta, TracksRun.Meta):
        unique_together = ("feature", "hash")

    @classmethod
    def get_or_create(cls, feature, value):
        # simple values: (int, float, str, bool, datetime)
        if not isinstance(value, dict):
            hash = hash_string(str(value))
        else:
            hash = hash_dict(value)
        try:
            return (
                cls.objects.create(feature=feature, value=value, hash=hash),
                False,
            )
        except IntegrityError:
            return cls.objects.get(feature=feature, hash=hash), True


def suggest_categorical_for_str_iterable(
    iterable: Iterable[str], key: str = None
) -> str:
    c = pd.Categorical(iterable)
    message = ""
    if len(c.categories) < len(c):
        if key != "":
            key_note = f" for feature {key}"
        else:
            key_note = ""
        message = f"You have few permissible values{key_note}, consider dtype 'cat' instead of 'str'"
    return message


def categoricals_from_df(df: pd.DataFrame) -> dict:
    """Returns categorical columns."""
    string_cols = [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {
        col: df[col]
        for col in df.columns
        if isinstance(df[col].dtype, CategoricalDtype)
    }
    for key in string_cols:
        message = suggest_categorical_for_str_iterable(df[key], key)
        if message:
            logger.warning(message)
    return categoricals
