import lnschema_bionty as bionty
from lamin_logger import logger
from typing_extensions import Literal

from lamindb.dev.db._add import add
from lamindb.dev.db._select import select
from lamindb.schema._table import Table


def add_features_and_featureset(
    features: dict,  # what is a features dict? can we have something more typed?
    feature_entity: Literal["gene", "protein", "cell_marker"],
    species: bionty.species,
    name: str = None,
):
    """Add features and featureset given a dictionary of input features.

    Test whether featureset of same name exists.

    If not, add all features that are in features_dict but not yet in the db.

    Link all features in the provided features_dict to the new featureset.
    """
    # check whether featureset already exists
    if name is not None:
        featureset = select(
            bionty.featureset, feature_entity=feature_entity, name=name
        ).one_or_none()
        if featureset is not None:
            logger.warning(f"featureset {name} already exists!")
            return featureset

    # get the id field of feature entity
    feature_id = features[next(iter(features))].keys()[-1]  # quite hard to read
    model = Table.get_model(feature_entity)
    exist_features = select(model, species_id=species).all()
    # only ingest the new features but link all features to the featureset
    exist_feature_keys, exist_feature_ids = set(), set()
    for (
        feature
    ) in exist_features:  # why do we have to write a for-loop here? rather pandas?
        exist_feature_keys.add(feature.__getattribute__(feature_id))
        exist_feature_ids.add(feature.id)
    # add non-exist features to the feature table
    new_features = []
    for k, v in features.items():
        if k in exist_feature_keys:
            continue
        new_features.append(model(**v))
    added = add(new_features)
    if not isinstance(added, list):
        added = list(added)
    feature_ids = [row.id for row in added] + list(exist_feature_ids)

    # add a featureset to the featureset table
    featureset = add(bionty.featureset(feature_entity=feature_entity, name=name))

    for feature_id in feature_ids:
        kwargs = {
            "featureset_id": featureset.id,
            f"{feature_entity}_id": feature_id,
        }
        add(getattr(bionty, f"featureset_{feature_entity}")(**kwargs))

    return featureset
