from __future__ import annotations

from typing import TYPE_CHECKING, Iterable, Type

import lamindb_setup as ln_setup
import numpy as np
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core.hashing import hash_set
from lnschema_core import Feature, FeatureSet, Record, ids
from lnschema_core.types import FieldAttr, ListLike

from lamindb._utils import attach_func_to_class_method

from ._feature import convert_numpy_dtype_to_lamin_feature_type
from ._record import init_self_from_db
from .core.exceptions import ValidationError
from .core.schema import (
    dict_related_model_to_related_name,
    get_related_name,
)

if TYPE_CHECKING:
    import pandas as pd

    from ._query_set import QuerySet

NUMBER_TYPE = "number"
DICT_KEYS_TYPE = type({}.keys())  # type: ignore


def validate_features(features: list[Record]) -> Record:
    """Validate and return feature type."""
    try:
        if len(features) == 0:
            raise ValueError("Provide list of features with at least one element")
    except TypeError:
        raise ValueError(
            "Please pass a ListLike of features, not a single feature"
        ) from None
    if not hasattr(features, "__getitem__"):
        raise TypeError("features has to be list-like")
    if not isinstance(features[0], Record):
        raise TypeError(
            "features has to store feature records! use .from_values() otherwise"
        )
    feature_types = {feature.__class__ for feature in features}
    if len(feature_types) > 1:
        raise TypeError("feature_set can only contain a single type")
    for feature in features:
        if feature._state.adding:
            raise ValueError("Can only construct feature sets from validated features")
    return next(iter(feature_types))  # return value in set of cardinality 1


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(FeatureSet, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 1:
        raise ValueError("Only one non-keyword arg allowed: features")
    features: Iterable[Record] = kwargs.pop("features") if len(args) == 0 else args[0]
    dtype: str | None = kwargs.pop("dtype") if "dtype" in kwargs else None
    name: str | None = kwargs.pop("name") if "name" in kwargs else None
    if len(kwargs) > 0:
        raise ValueError("Only features, dtype, name are valid keyword arguments")
    # now code
    features_registry = validate_features(features)
    if dtype is None:
        dtype = None if features_registry == Feature else NUMBER_TYPE
    n_features = len(features)
    features_hash = hash_set({feature.uid for feature in features})
    feature_set = FeatureSet.filter(hash=features_hash).one_or_none()
    if feature_set is not None:
        logger.debug(f"loaded: {feature_set}")
        init_self_from_db(self, feature_set)
        return None
    else:
        hash = features_hash
    self._features = (get_related_name(features_registry), features)

    super(FeatureSet, self).__init__(
        uid=ids.base62_20(),
        name=name,
        dtype=get_type_str(dtype),
        n=n_features,
        registry=features_registry.__get_name_with_schema__(),
        hash=hash,
    )


@doc_args(FeatureSet.save.__doc__)
def save(self, *args, **kwargs) -> FeatureSet:
    """{}"""  # noqa: D415
    super(FeatureSet, self).save(*args, **kwargs)
    if hasattr(self, "_features"):
        related_name, records = self._features
        getattr(self, related_name).set(records)
    return self


def get_type_str(dtype: str | None) -> str | None:
    if dtype is not None:
        type_str = dtype.__name__ if not isinstance(dtype, str) else dtype  # type: ignore
    else:
        type_str = None
    return type_str


@classmethod  # type:ignore
@doc_args(FeatureSet.from_values.__doc__)
def from_values(
    cls,
    values: ListLike,
    field: FieldAttr = Feature.name,
    type: str | None = None,
    name: str | None = None,
    mute: bool = False,
    organism: Record | str | None = None,
    source: Record | None = None,
    raise_validation_error: bool = True,
) -> FeatureSet:
    """{}"""  # noqa: D415
    if not isinstance(field, FieldAttr):
        raise TypeError("Argument `field` must be a Record field, e.g., `Feature.name`")
    if len(values) == 0:
        raise ValueError("Provide a list of at least one value")
    if isinstance(values, DICT_KEYS_TYPE):
        values = list(values)
    registry = field.field.model
    if registry != Feature and type is None:
        type = NUMBER_TYPE
        logger.debug("setting feature set to 'number'")
    validated = registry.validate(values, field=field, mute=mute, organism=organism)
    values_array = np.array(values)
    validated_values = values_array[validated]
    if validated.sum() != len(values):
        not_validated_values = values_array[~validated]
        msg = (
            f"These values could not be validated: {not_validated_values.tolist()}\n"
            f"If there are no typos, add them to their registry: {registry}"
        )
        if raise_validation_error:
            raise ValidationError(msg)
        elif len(validated_values) == 0:
            return None  # temporarily return None here
    validated_features = registry.from_values(
        validated_values,
        field=field,
        organism=organism,
        source=source,
    )
    feature_set = FeatureSet(
        features=validated_features,
        name=name,
        dtype=get_type_str(type),
    )
    return feature_set


@classmethod  # type:ignore
@doc_args(FeatureSet.from_df.__doc__)
def from_df(
    cls,
    df: pd.DataFrame,
    field: FieldAttr = Feature.name,
    name: str | None = None,
    mute: bool = False,
    organism: Record | str | None = None,
    source: Record | None = None,
) -> FeatureSet | None:
    """{}"""  # noqa: D415
    registry = field.field.model
    validated = registry.validate(df.columns, field=field, mute=mute, organism=organism)
    if validated.sum() == 0:
        if mute is True:
            logger.warning("no validated features, skip creating feature set")
        return None
    if registry == Feature:
        validated_features = Feature.from_df(df.loc[:, validated])
        feature_set = FeatureSet(validated_features, name=name, dtype=None)
    else:
        dtypes = [col.dtype for (_, col) in df.loc[:, validated].items()]
        if len(set(dtypes)) != 1:
            raise ValueError(f"data types are heterogeneous: {set(dtypes)}")
        dtype = convert_numpy_dtype_to_lamin_feature_type(dtypes[0])
        validated_features = registry.from_values(
            df.columns[validated],
            field=field,
            organism=organism,
            source=source,
        )
        feature_set = FeatureSet(
            features=validated_features,
            name=name,
            dtype=get_type_str(dtype),
        )
    return feature_set


@property  # type: ignore
@doc_args(FeatureSet.members.__doc__)
def members(self) -> QuerySet:
    """{}"""  # noqa: D415
    if self._state.adding:
        # this should return a queryset and not a list...
        # need to fix this
        return self._features[1]
    related_name = self._get_related_name()
    if related_name is None:
        related_name = "features"
    return self.__getattribute__(related_name).all()


def _get_related_name(self: FeatureSet) -> str:
    feature_sets_related_models = dict_related_model_to_related_name(self)
    related_name = feature_sets_related_models.get(self.registry)
    return related_name


METHOD_NAMES = [
    "__init__",
    "from_values",
    "from_df",
    "save",
]

if ln_setup._TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(FeatureSet, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, FeatureSet, globals())

FeatureSet.members = members
FeatureSet._get_related_name = _get_related_name
