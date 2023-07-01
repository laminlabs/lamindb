from typing import List, Optional

from lamin_logger import logger
from lnschema_core import ORM, Feature, FeatureSet

from lamindb.dev.hashing import hash_set

from ._from_values import Field, ListLike, get_or_create_records, index_iterable


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
    if not isinstance(features, ListLike):
        raise TypeError("features has to be list-like")
    if not isinstance(features[0], ORM):
        raise TypeError(
            "features has to store feature records! use .from_values() otherwise"
        )
    feature_types = set([feature.model for feature in features])
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
    assert len(kwargs) == 0
    features_type = validate_features(features)
    related_name = get_related_name(features_type)
    features_hash = hash_set({feature.id for feature in features})
    self._features = (related_name, features)
    field_name = "id"
    super(FeatureSet, self).__init__(
        id=features_hash, type=features_type.__name__, field=field_name
    )


def save(self, *args, **kwargs):
    super(FeatureSet, self).save(*args, **kwargs)
    if hasattr(self, "_features"):
        related_name, records = self._features
        getattr(self, related_name).set(records)


@classmethod  # type:ignore
def from_values(
    cls, values: ListLike, field: Field = Feature.name, species: Optional[str] = None
) -> FeatureSet:
    if not isinstance(field, Field):
        raise TypeError("Argument `field` must be an ORM field, e.g., `Feature.name`")
    orm = field.field.model
    iterable_idx = index_iterable(values)
    features_hash = hash_set(set(iterable_idx))
    feature_set = FeatureSet.select(id=features_hash, type=orm.__name__).one_or_none()
    if feature_set is not None:
        logger.info("Returning an existing feature_set")
    else:
        from_bionty = orm.__module__.startswith("lnschema_bionty.")
        records = get_or_create_records(
            iterable=iterable_idx,
            field=field,
            from_bionty=from_bionty,
            species=species,
        )
        feature_set = FeatureSet(
            id=features_hash,
            type=orm.__name__,
            field=field.field_name,
            features=records,
        )
    return feature_set


FeatureSet.__init__ = __init__
FeatureSet.save = save
FeatureSet.from_values = from_values
