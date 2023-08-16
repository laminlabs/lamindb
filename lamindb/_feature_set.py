from typing import Iterable, List, Optional, Type, Union

import pandas as pd
from lamin_utils import logger
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Feature, FeatureSet, Modality, Registry, ids
from lnschema_core.types import FieldAttr, ListLike

from lamindb.dev.hashing import hash_set
from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import get_or_create_records, index_iterable
from ._registry import init_self_from_db
from ._save import bulk_create


def get_related_name(features_type: Registry):
    candidates = [
        field.related_name
        for field in FeatureSet._meta.related_objects
        if field.related_model == features_type
    ]
    if not candidates:
        raise ValueError(
            f"Can't create feature sets from {features_type.__name__} because it's not"
            " related to it!\nYou need to create a link model between FeatureSet and"
            " your Registry in your custom schema.\nTo do so, add a"
            " line:\nfeature_sets = models.ManyToMany(FeatureSet,"
            " related_name='mythings')\n"
        )
    return candidates[0]


def sanity_check_features(features: List[Registry]) -> Registry:
    """Validate and return feature type."""
    if len(features) == 0:
        raise ValueError("provide list of features with at least one element")
    if not hasattr(features, "__getitem__"):
        raise TypeError("features has to be list-like")
    if not isinstance(features[0], Registry):
        raise TypeError(
            "features has to store feature records! use .from_values() otherwise"
        )
    feature_types = set([feature.__class__ for feature in features])
    if len(feature_types) > 1:
        raise ValueError("feature_set can only contain a single type")
    return next(iter(feature_types))  # return value in set of cardinality 1


def get_validated_features(
    features: List[Registry], field: FieldAttr
) -> List[Registry]:
    validated_features = []
    non_validated_features = []
    for feature in features:
        if feature._state.adding and not (
            hasattr(feature, "_from_bionty") and feature._from_bionty
        ):
            non_validated_features.append(getattr(feature, field.field.name))
        else:
            validated_features.append(feature)
    if non_validated_features:
        non_validated_features_display = ",".join(non_validated_features)
        logger.warning(
            f"ignoring non-validated features: {non_validated_features_display}"
        )
    return validated_features


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(FeatureSet, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 1:
        raise ValueError("Only one non-keyword arg allowed: features")
    features: Iterable[Registry] = kwargs.pop("features") if len(args) == 0 else args[0]
    type: Optional[Union[type, str]] = kwargs.pop("type") if "type" in kwargs else None
    modality: Optional[str] = kwargs.pop("modality") if "modality" in kwargs else None
    name: Optional[str] = kwargs.pop("name") if "name" in kwargs else None
    # hash is only internally used
    hash: Optional[str] = kwargs.pop("hash") if "hash" in kwargs else None
    if len(kwargs) > 0:
        raise ValueError(
            "Only features, type, modality, name are valid keyword arguments"
        )

    # now code
    features_orm = sanity_check_features(features)
    if features_orm == Feature:
        type = None
    else:
        type = float
    n_features = len(features)
    if hash is None:
        features_hash = hash_set({feature.id for feature in features})
        feature_set = FeatureSet.filter(hash=features_hash).one_or_none()
        if feature_set is not None:
            logger.success(f"loaded: {feature_set}")
            init_self_from_db(self, feature_set)
            return None
        else:
            hash = features_hash
    self._features = (get_related_name(features_orm), features)
    if type is not None:
        type_str = type.__name__ if not isinstance(type, str) else type
    else:
        type_str = None
    if modality is not None:
        if isinstance(modality, str):
            modality_record = Modality.filter(name=modality).one_or_none()
            if modality_record is None:
                modality_record = Modality(name=modality)
                modality_record.save()
        elif isinstance(modality, Modality):
            modality_record = modality
        else:
            raise ValueError("modality needs to be string or Modality record")
    else:
        modality_record = modality
    super(FeatureSet, self).__init__(
        id=ids.base62_20(),
        name=name,
        type=type_str,
        n=n_features,
        modality=modality_record,
        registry=features_orm.__get_name_with_schema__(),
        hash=hash,
    )


@doc_args(FeatureSet.save.__doc__)
def save(self, *args, **kwargs) -> None:
    """{}"""
    super(FeatureSet, self).save(*args, **kwargs)
    if hasattr(self, "_features"):
        related_name, records = self._features
        # if values are stored in their dedicated table, we can bulk_create
        if related_name != "features":
            bulk_create(records)
        else:
            # otherwise, we currently need to save one by one
            for record in records:
                record.save()
        getattr(self, related_name).set(records)


@classmethod  # type:ignore
@doc_args(FeatureSet.from_values.__doc__)
def from_values(
    cls,
    values: ListLike,
    field: FieldAttr = Feature.name,
    type: Optional[Union[Type, str]] = None,
    name: Optional[str] = None,
    modality: Optional[str] = None,
    **kwargs,
) -> Optional["FeatureSet"]:
    """{}"""
    if not isinstance(field, FieldAttr):
        raise TypeError(
            "Argument `field` must be a Registry field, e.g., `Feature.name`"
        )
    if len(values) == 0:
        raise ValueError("Provide a list of at least one value")
    registry = field.field.model
    if registry.__name__ == "Feature":
        raise ValueError("Please use from_df() instead of from_values()")
    iterable_idx = index_iterable(values)
    if not isinstance(iterable_idx[0], (str, int)):
        raise TypeError("values should be list-like of str or int")
    from_bionty = registry.__module__.startswith("lnschema_bionty")
    features = get_or_create_records(
        iterable=iterable_idx,
        field=field,
        from_bionty=from_bionty,
        **kwargs,
    )
    validated_features = get_validated_features(features, field)
    validated_feature_ids = [feature.id for feature in validated_features]
    features_hash = hash_set(set(validated_feature_ids))
    feature_set = FeatureSet.filter(hash=features_hash).one_or_none()
    if feature_set is not None:
        logger.success(f"loaded {feature_set}")
    else:
        if type is not None:
            type_str = type.__name__ if not isinstance(type, str) else type
        else:
            type_str = None
        if validated_features:
            feature_set = FeatureSet(
                features=validated_features,
                hash=features_hash,
                name=name,
                modality=modality,
                type=type_str,
            )
        else:
            feature_set = None
    return feature_set


@classmethod  # type:ignore
@doc_args(FeatureSet.from_df.__doc__)
def from_df(
    cls,
    df: "pd.DataFrame",
    name: Optional[str] = None,
) -> Optional["FeatureSet"]:
    """{}"""
    features = Feature.from_df(df)
    validated_features = get_validated_features(features, Feature.name)
    if validated_features:
        feature_set = FeatureSet(validated_features, name=name)
    else:
        logger.warning("no validated features, skip creating feature set")
        feature_set = None
        # raise ValidationError("Dataframe columns contain no validated feature names")
    return feature_set


METHOD_NAMES = [
    "__init__",
    "from_values",
    "from_df",
    "save",
]

if _TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(FeatureSet, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, FeatureSet, globals())
