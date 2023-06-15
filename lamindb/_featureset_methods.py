from typing import Dict, Iterable, Optional

from django.db import models
from lamin_logger import logger
from lnschema_core import FeatureSet

from lamindb._select import select
from lamindb.dev.hashing import hash_set

from ._parse import Field, ListLike, get_or_create_records, index_iterable


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

    featureset = select(
        FeatureSet,
        id=features_hash,
        type=related_name,
    ).one_or_none()
    if featureset is not None:
        logger.info("Returning an existing featureset")
    else:
        records = get_or_create_records(
            iterable=iterable_idx, field=field, species=species
        )
        featureset = FeatureSet(
            id=features_hash, type=related_name, **{related_name: records}
        )
    return featureset


@classmethod  # type: ignore
def from_iterable(
    cls,
    iterable: Iterable,
    field: models.CharField,
    species: str = None,
):
    """Parse iterable & return featureset & records."""
    from lamindb._featureset_methods import parse_features_from_iterable

    features = parse_features_from_iterable(
        iterable=iterable,
        field=field,
        species=species,
    )
    return features


def __init__(featureset, *args, **kwargs):  # type: ignore
    related_names = [i.related_name for i in featureset.__class__._meta.related_objects]

    relationships: Dict = {}
    for related_name in related_names:
        if related_name in kwargs:
            relationships[related_name] = kwargs.pop(related_name)
    featureset._relationships = relationships

    super(FeatureSet, featureset).__init__(*args, **kwargs)


def save(featureset, *args, **kwargs):
    super(FeatureSet, featureset).save(*args, **kwargs)
    for key, records in featureset._relationships.items():
        [r.save() for r in records]
        getattr(featureset, key).set(records)


FeatureSet.from_iterable = from_iterable
FeatureSet.__init__ = __init__
FeatureSet.save = save
