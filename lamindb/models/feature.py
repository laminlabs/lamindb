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
from lamindb.errors import FieldValidationError, ValidationError

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


def parse_dtype(dtype_str: str, is_param: bool = False) -> list[dict[str, str]]:
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
    split_result = dtype_str.split("[")
    # has sub type
    sub_type_str = ""
    if len(split_result) == 2:
        registry_str = split_result[0]
        assert "]" in split_result[1]  # noqa: S101
        sub_type_field_split = split_result[1].split("].")
        if len(sub_type_field_split) == 1:
            sub_type_str = sub_type_field_split[0].strip("]")
            field_str = ""
        else:
            sub_type_str = sub_type_field_split[0]
            field_str = sub_type_field_split[1]
    elif len(split_result) == 1:
        registry_field_split = split_result[0].split(".")
        if (
            len(registry_field_split) == 2 and registry_field_split[1][0].isupper()
        ) or len(registry_field_split) == 3:
            # bionty.CellType or bionty.CellType.name
            registry_str = f"{registry_field_split[0]}.{registry_field_split[1]}"
            field_str = (
                "" if len(registry_field_split) == 2 else registry_field_split[2]
            )
        else:
            # ULabel or ULabel.name
            registry_str = registry_field_split[0]
            field_str = (
                "" if len(registry_field_split) == 1 else registry_field_split[1]
            )
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
                raise ImportError(
                    f"Can not parse dtype {dtype_str} because {module_name_attempt} "
                    f"was not found.\nInstall the module with `pip install {module_name_attempt}`\n"
                    "and also add the module to this instance via instance settings page "
                    "under 'schema modules'."
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
    return {
        "registry": registry,  # should be typed as CanCurate
        "registry_str": registry_str,
        "subtype_str": sub_type_str,
        "field_str": field_str,
        "field": getattr(registry, field_str),
    }


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


def process_init_feature_param(args, kwargs, is_param: bool = False):
    # now we proceed with the user-facing constructor
    if len(args) != 0:
        raise ValueError("Only keyword args allowed")
    name: str = kwargs.pop("name", None)
    dtype: type | str | None = kwargs.pop("dtype", None)
    is_type: bool = kwargs.pop("is_type", None)
    type_: Feature | str | None = kwargs.pop("type", None)
    description: str | None = kwargs.pop("description", None)
    if kwargs:
        valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Feature)])
        raise FieldValidationError(f"Only {valid_keywords} are valid keyword arguments")
    kwargs["name"] = name
    kwargs["type"] = type_
    kwargs["is_type"] = is_type
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

    A feature often represents a dimension of a dataset, such as a column in a
    `DataFrame`. The `Feature` registry organizes metadata of features.

    The `Feature` registry helps you organize and query datasets based on their
    features and corresponding label annotations. For instance, when working
    with a "T cell" label, it could be measured through different features
    such as `"cell_type_by_expert"` where an expert manually classified the
    cell, or `"cell_type_by_model"` where a computational model made the
    classification.

    The two most important metadata of a feature are its `name` and the `dtype`.
    In addition to typical data types, LaminDB has a `"num"` `dtype` to
    concisely denote the union of all numerical types.

    Args:
        name: `str` Name of the feature, typically.  column name.
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

        A simple `"str"` feature.

        >>> ln.Feature(
        ...     name="sample_note",
        ...     dtype="str",
        ... ).save()

        A dtype `"cat[ULabel]"` can be more easily passed as below.

        >>> ln.Feature(
        ...     name="project",
        ...     dtype=ln.ULabel,
        ... ).save()

        A dtype `"cat[ULabel|bionty.CellType]"` can be more easily passed as below.

        >>> ln.Feature(
        ...     name="cell_type",
        ...     dtype=[ln.ULabel, bt.CellType],
        ... ).save()

        A multivalue feature with a list of cell types.

        >>> ln.Feature(
        ...     name="cell_types",
        ...     dtype=list[bt.CellType],  # or list[str] for a list of strings
        ... ).save()

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
        dtype = kwargs.get("dtype", None)
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
            assert "|" not in dtype_str  # noqa: S101
            assert "]]" not in dtype_str  # noqa: S101
            fill_in = ", ".join(
                f"{key}='{value}'" for (key, value) in cat_filters.items()
            )
            dtype_str = dtype_str.replace("]", f"[{fill_in}]]")
            self.dtype = dtype_str
        if not self._state.adding:
            if not (
                self.dtype.startswith("cat")
                if dtype == "cat"
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
