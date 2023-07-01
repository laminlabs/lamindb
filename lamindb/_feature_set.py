from typing import Dict, Optional, Union

from lamin_logger import logger
from lnschema_core import Feature, FeatureSet

from lamindb._select import select
from lamindb.dev.hashing import hash_set

from ._from_values import Field, ListLike, get_or_create_records, index_iterable


# expose to user via ln.FeatureSet
def parse_features_from_iterable(
    iterable: ListLike,
    field: Field,
    species: Optional[str] = None,
):
    # get related_name of the field class from FeatureSet class
    model = field.field.model
    related_name = [
        i.related_name
        for i in FeatureSet._meta.related_objects
        if i.related_model == model
    ]
    if len(related_name) == 0:
        raise AssertionError(
            f"Can't create featuresets from {model.__name__}! Check your schema!"
        )
    else:
        related_name = related_name[0]

    iterable_idx = index_iterable(iterable)

    features_hash = hash_set(set(iterable_idx))

    feature_set = select(
        FeatureSet,
        id=features_hash,
        type=related_name,
    ).one_or_none()

    from_bionty = (
        True if field.field.model.__module__.startswith("lnschema_bionty.") else False
    )
    if feature_set is not None:
        logger.info("Returning an existing feature_set")
    else:
        if species is not None:
            kwargs = dict(species=species)
        else:
            kwargs = dict()
        records = get_or_create_records(
            iterable=iterable_idx, field=field, from_bionty=from_bionty, **kwargs
        )
        feature_set = FeatureSet(
            id=features_hash, type=related_name, **{related_name: records}
        )
    return feature_set


def __init__(self, *args, **kwargs):  # type: ignore
    related_names = [i.related_name for i in self.__class__._meta.related_objects]

    relationships: Dict = {}
    for related_name in related_names:
        if related_name in kwargs:
            relationships[related_name] = kwargs.pop(related_name)
    self._relationships = relationships

    super(FeatureSet, self).__init__(*args, **kwargs)


def save(feature_set, *args, **kwargs):
    super(FeatureSet, feature_set).save(*args, **kwargs)
    for key, records in feature_set._relationships.items():
        [r.save() for r in records]
        getattr(feature_set, key).set(records)


@classmethod  # type:ignore
def from_values(
    cls, values: ListLike, field: Union[Field, str] = Feature.name, **kwargs
):
    if isinstance(field, str):
        field = getattr(cls, field)
    if not isinstance(field, Field):  # field is DeferredAttribute
        raise TypeError(
            "field must be a string or an ORM field, e.g., `CellType.name`!"
        )
    feature_set = parse_features_from_iterable(
        iterable=values,
        field=field,
        species=kwargs.get("species"),
    )
    return feature_set


FeatureSet.__init__ = __init__
FeatureSet.save = save
FeatureSet.from_values = from_values
