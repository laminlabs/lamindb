from typing import Any

import pandas as pd
from lnschema_core import Features

from lamindb._select import select
from lamindb.dev.hashing import hash_set

from ._parse import Field, ListLike, get_or_create_records


def parse_features(df: pd.DataFrame, bionty_object: Any, **map_kwargs) -> None:
    """Link features to a knowledge table.

    Args:
        df: a DataFrame
        bionty_object: Features reference class, bionty.{entity}()
    """
    from bionty import CellMarker, Gene, Protein

    if "__field__" in map_kwargs:
        field = map_kwargs.pop("__field__")
        map_kwargs["reference_id"] = field.name
    df_curated = bionty_object.curate(df=df, **map_kwargs)
    # ._parsing_id only exist and set to be the reference_id after `.curate`` is called
    parsing_id = bionty_object._parsing_id

    # logging of curation
    n = df_curated["__curated__"].count()
    n_mapped = df_curated["__curated__"].sum()
    log = {  # noqa  TODO: store this somewhere in the db
        "feature": parsing_id,
        "n_mapped": n_mapped,
        "percent_mapped": round(n_mapped / n * 100, 1),
        "unmapped": df_curated.index[~df_curated["__curated__"]],
    }

    features_hash = hash_set(set(df_curated.index))

    features_type = bionty_object._entity

    features = select(
        Features,
        id=features_hash,
        type=features_type,
    ).one_or_none()
    if features is not None:
        return features  # features already exists!
    else:
        features = Features(id=features_hash, type=features_type)
        records = get_or_create_records(
            iterable=df_curated.index,
            field=field,
            species=bionty_object.species,
        )

        if isinstance(bionty_object, Gene):
            for record in records:
                features.genes.append(record)
        elif isinstance(bionty_object, Protein):
            for record in records:
                features.proteins.append(record)
        elif isinstance(bionty_object, CellMarker):
            for record in records:
                features.cell_markers.append(record)
        else:
            raise NotImplementedError(
                "Features parsing is only supported for 'Gene', 'Protein' and"
                " 'CellMarker'!"
            )

    return features


# expose to user via ln.Features
def parse_features_from_iterable(
    iterable: ListLike,
    field: Field,
    species: str = None,
):
    # No entries are made for NAs, '', None
    iterable = [i for i in set(iterable) if not (pd.isnull(i) or i == "" or i == " ")]
    entity = field.field.model
    reference = entity.bionty(species=species)
    df = pd.DataFrame(index=list(iterable))
    return parse_features(df, reference)
