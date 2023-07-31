from typing import Iterable, List, Optional, Type, Union

import pandas as pd
from django.db.models.query_utils import DeferredAttribute as Field
from lamin_utils import logger
from lamindb_setup.dev._docs import doc_args
from lnschema_core import ORM, Feature, FeatureSet, Modality, ids
from lnschema_core.types import ListLike

from lamindb.dev.hashing import hash_set
from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import get_or_create_records, index_iterable
from ._orm import init_self_from_db
from ._save import bulk_create


def get_related_name(features_type: ORM):
    candidates = [
        field.related_name
        for field in FeatureSet._meta.related_objects
        if field.related_model == features_type
    ]
    if not candidates:
        raise ValueError(
            f"Can't create feature sets from {features_type.__name__} because it's not"
            " related to it!\nYou need to create a link model between FeatureSet and"
            " your ORM in your custom schema.\nTo do so, add a line:\nfeature_sets ="
            " models.ManyToMany(FeatureSet, related_name='mythings')\n"
        )
    return candidates[0]


def validate_features(features: List[ORM]) -> ORM:
    """Validate and return feature type."""
    if len(features) == 0:
        raise ValueError("provide list of features with at least one element")
    if not hasattr(features, "__getitem__"):
        raise TypeError("features has to be list-like")
    if not isinstance(features[0], ORM):
        raise TypeError(
            "features has to store feature records! use .from_values() otherwise"
        )
    feature_types = set([feature.__class__ for feature in features])
    if len(feature_types) > 1:
        raise ValueError("feature_set can only contain a single type")
    return next(iter(feature_types))  # return value in set of cardinality 1


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(FeatureSet, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 1:
        raise ValueError("Only one non-keyword arg allowed: features")
    features: Iterable[ORM] = kwargs.pop("features") if len(args) == 0 else args[0]
    ref_field: Optional[str] = (
        kwargs.pop("ref_field") if "ref_field" in kwargs else "id"
    )
    type: Optional[Union[type, str]] = kwargs.pop("type") if "type" in kwargs else None
    modality: Optional[str] = kwargs.pop("modality") if "modality" in kwargs else None
    name: Optional[str] = kwargs.pop("name") if "name" in kwargs else None
    # hash is only internally used
    hash: Optional[str] = kwargs.pop("hash") if "hash" in kwargs else None
    if len(kwargs) > 0:
        raise ValueError(
            "Only features, ref_field, type, modality, name are valid keyword arguments"
        )

    # now code
    features_orm = validate_features(features)
    if features_orm == Feature:
        type = None
    else:
        type = float
    n_features = len(features)
    if hash is None:
        features_hash = hash_set({feature.id for feature in features})
        feature_set = FeatureSet.filter(hash=features_hash).one_or_none()
        if feature_set is not None:
            logger.info(f"Loaded {feature_set}")
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
        ref_field=f"{features_orm.__get_name_with_schema__()}.{ref_field}",
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
    field: Field = Feature.name,
    type: Optional[Union[Type, str]] = None,
    name: Optional[str] = None,
    modality: Optional[str] = None,
    **kwargs,
) -> "FeatureSet":
    """{}"""
    if not isinstance(field, Field):
        raise TypeError("Argument `field` must be an ORM field, e.g., `Feature.name`")
    if len(values) == 0:
        raise ValueError("Provide a list of at least one value")
    ORM = field.field.model
    if isinstance(ORM, Feature):
        raise ValueError("Please use from_df() instead of from_values()")
    iterable_idx = index_iterable(values)
    if not isinstance(iterable_idx[0], (str, int)):
        raise TypeError("values should be list-like of str or int")
    features_hash = hash_set(set(iterable_idx))
    feature_set = FeatureSet.filter(hash=features_hash).one_or_none()
    if feature_set is not None:
        logger.info(f"Loaded {feature_set}")
    else:
        from_bionty = ORM.__module__.startswith("lnschema_bionty")
        records = get_or_create_records(
            iterable=iterable_idx,
            field=field,
            from_bionty=from_bionty,
            **kwargs,
        )
        # type_str = type.__name__ if not isinstance(type, str) else type
        feature_set = FeatureSet(
            features=records,
            hash=features_hash,
            name=name,
            modality=modality,
            type=type,
            ref_field=field.field.name,
        )
    return feature_set


@classmethod  # type:ignore
@doc_args(FeatureSet.from_df.__doc__)
def from_df(
    cls,
    df: "pd.DataFrame",
    name: Optional[str] = None,
) -> "FeatureSet":
    """{}"""
    features = Feature.from_df(df)
    feature_set = FeatureSet(features, name=name)
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
