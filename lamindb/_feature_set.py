from typing import List, Optional

from django.db.models.query_utils import DeferredAttribute as Field
from lamin_logger import logger
from lamindb_setup.dev._docs import doc_args
from lnschema_core import ORM, Feature, FeatureSet
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
    features: List[ORM] = kwargs.pop("features") if len(args) == 0 else args[0]
    field: Optional[str] = kwargs.pop("field") if "field" in kwargs else None
    id: Optional[str] = kwargs.pop("id") if "id" in kwargs else None
    features_type = validate_features(features)
    related_name = get_related_name(features_type)
    if id is None:
        features_hash = hash_set({feature.id for feature in features})
        feature_set = FeatureSet.select(id=features_hash).one_or_none()
        if feature_set is not None:
            logger.info("Returning an existing feature_set")
            init_self_from_db(self, feature_set)
            return None
        else:
            id = features_hash
    self._features = (related_name, features)
    if field is None:
        field = "id"
    super(FeatureSet, self).__init__(
        id=id, type=features_type.__name_with_type__(), field=field
    )


@doc_args(FeatureSet.save.__doc__)
def save(self, *args, **kwargs) -> None:
    """{}"""
    super(FeatureSet, self).save(*args, **kwargs)
    if hasattr(self, "_features"):
        related_name, records = self._features
        bulk_create(records)
        getattr(self, related_name).set(records)


@classmethod  # type:ignore
@doc_args(FeatureSet.from_values.__doc__)
def from_values(
    cls, values: ListLike, field: Field = Feature.name, **kwargs
) -> "FeatureSet":
    """{}"""
    if not isinstance(field, Field):
        raise TypeError("Argument `field` must be an ORM field, e.g., `Feature.name`")
    if len(values) == 0:
        raise ValueError("Provide a list of at least one value")
    orm = field.field.model
    iterable_idx = index_iterable(values)
    if not isinstance(iterable_idx[0], (str, int)):
        raise TypeError("values should be list-like of str or int")
    features_hash = hash_set(set(iterable_idx))
    feature_set = FeatureSet.select(id=features_hash).one_or_none()
    if feature_set is not None:
        logger.info("Returning an existing feature_set")
    else:
        from_bionty = orm.__module__.startswith("lnschema_bionty")
        records = get_or_create_records(
            iterable=iterable_idx,
            field=field,
            from_bionty=from_bionty,
            **kwargs,
        )
        feature_set = FeatureSet(
            id=features_hash,
            field=field.field.name,
            features=records,
        )
    return feature_set


METHOD_NAMES = [
    "__init__",
    "from_values",
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
