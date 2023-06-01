from typing import Any

import pandas as pd
from lnschema_core import Features

from lamindb.dev.db._select import select
from lamindb.dev.hashing import hash_set
from lamindb.dev.storage import load_to_memory

from ._parse import get_or_create_records


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


def get_features(bionty_object, iterable=None, file_privates=None, **map_kwargs):
    """Updates file in place."""
    if file_privates is not None:
        memory_rep = file_privates["_memory_rep"]
        if memory_rep is None:
            memory_rep = load_to_memory(file_privates["_local_filepath"])
        try:
            df = getattr(memory_rep, "var")  # for AnnData objects
            if callable(df):
                df = memory_rep
        except AttributeError:
            df = memory_rep
    elif iterable is not None:
        df = pd.DataFrame(index=list(iterable))
    else:
        raise KeyError

    return parse_features(df, bionty_object, **map_kwargs)
