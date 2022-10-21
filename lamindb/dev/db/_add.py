from typing import List, Union

import lnschema_bionty as bionty
import sqlmodel as sqm
from lamin_logger import logger
from lndb_setup import settings
from typing_extensions import Literal

from lamindb.schema._table import Table

from ._select import select


def add(*rows: sqm.SQLModel) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
    """Insert or update rows.

    Inserts a new row if it doesn't exist. Updates the row if it exists already.

    To achieve the latter, query the `row` with `.get` or `.select` before
    passing it to add.

    Guide: :doc:`/db/guide/add-delete`.

    Example:

    >>> experiment = ln.get(wetlab.experiment, "93jIJFla")
    >>> experiment.name = "New name"
    >>> add(experiment)

    Args:
        *rows: One or multiple rows as instances of `SQLModel`.
    """
    with sqm.Session(settings.instance.db_engine()) as session:
        for row in rows:
            session.add(row)
        session.commit()
        for row in rows:
            session.refresh(row)
    settings.instance._update_cloud_sqlite_file()
    if len(rows) > 1:
        return rows
    else:
        return rows[0]


def add_featureset_from_features(
    features_dict: dict,
    feature_entity: Literal["gene", "protein", "cell_marker"],
    species: bionty.species,
    featureset_name: str = None,
):
    """Insert a featureset."""
    # check whether featureset already exists
    if featureset_name is not None:
        featureset_result = (
            select(bionty.featureset)
            .where(
                bionty.featureset.feature_entity == feature_entity,
                bionty.featureset.name == featureset_name,
            )
            .one_or_none()
        )
        if featureset_result is not None:
            logger.warning(f"Featureset {featureset_name} already exists!")
            return featureset_result

    # get the id field of feature entity
    feature_id = features_dict[next(iter(features_dict))].keys()[-1]
    model = Table.get_model(feature_entity)
    allfeatures = select(model).where(model.species_id == species.id).all()  # type: ignore  # noqa
    # only ingest the new features but link all features to the featureset
    exist_feature_keys = set()
    exist_feature_ids = set()
    for feature in allfeatures:
        exist_feature_keys.add(feature.__getattribute__(feature_id))
        exist_feature_ids.add(feature.id)

    # add a featureset to the featureset table
    featureset = add.featureset(  # type: ignore
        feature_entity=feature_entity, name=featureset_name
    )

    # add features to the feature table
    kwargs_list = []
    for k, v in features_dict.items():
        if k in exist_feature_keys:
            continue
        kwargs_list.append(v)
    added = add(kwargs_list)
    if not isinstance(added, list):
        added = list(added)
    feature_ids = added + list(exist_feature_ids)
    for feature_id in feature_ids:
        kwargs = {
            "featureset_id": featureset.id,
            f"{feature_entity}_id": feature_id,
        }
        _ = getattr(add, f"featureset_{feature_entity}")(**kwargs)

    return featureset
