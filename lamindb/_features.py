from typing import Optional

from lamin_logger import logger
from lnschema_core import Featureset

from lamindb._select import select
from lamindb.dev.hashing import hash_set

from ._parse import Field, ListLike, get_or_create_records, index_iterable


# expose to user via ln.Featureset
def parse_features_from_iterable(
    iterable: ListLike,
    field: Field,
    species: Optional[str] = None,
):
    # get related_name of the field class from Featureset class
    model = field.field.model
    related_name = [
        i.related_name
        for i in Featureset._meta.related_objects
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
        Featureset,
        id=features_hash,
        type=related_name,
    ).one_or_none()
    if featureset is not None:
        logger.info("Returning an existing featureset")
    else:
        records = get_or_create_records(
            iterable=iterable_idx, field=field, species=species
        )
        featureset = Featureset(
            id=features_hash, type=related_name, **{related_name: records}
        )
    return featureset
