from typing import Any, List

import pandas as pd
from lndb_storage import load_to_memory
from lnschema_core import Features

from lamindb.dev.db._add import add
from lamindb.dev.db._select import select
from lamindb.dev.hashing import hash_set
from lamindb.schema._table import table_meta


def get_features_records(
    parsing_id: str,
    features_ref: Any,
    df_curated: pd.DataFrame,
) -> List:
    import lnschema_bionty as bionty

    # insert species entry if not exists
    species = select(bionty.Species, name=features_ref.species).one_or_none()
    if species is None:
        species = add(bionty.Species(name=features_ref.species))

    model = table_meta.get_model(f"bionty.{features_ref.entity}")

    # all existing feature records of the species in the db
    stmt = (
        select(model)
        .where(getattr(model, parsing_id).in_(df_curated.index))
        .where(getattr(model, "species_id") == species.id)  # type:ignore
    )
    records = stmt.all()
    records_df = df_curated.index.intersection(stmt.df()[parsing_id])

    # new records to be appended
    new_ids = df_curated.index.difference(records_df)
    if len(new_ids) > 0:
        # mapped new_ids
        reference_df = features_ref.df.set_index(features_ref.reference_id)
        mapped = reference_df.loc[reference_df.index.intersection(new_ids)].copy()
        mapped.index.name = parsing_id
        if mapped.shape[0] > 0:
            for kwargs in mapped.reset_index().to_dict(orient="records"):
                kwargs["species_id"] = species.id  # type:ignore
                record = model(**kwargs)
                records.append(record)
        # unmapped new_ids
        unmapped = set(new_ids).difference(mapped.index)
        if len(unmapped) > 0:
            for i in unmapped:
                record = model(
                    **{parsing_id: i, "species_id": species.id}  # type:ignore
                )
                records.append(record)
    return records


def parse_features(df: pd.DataFrame, features_ref: Any, **curate_kwargs) -> None:
    """Link features to a knowledge table.

    Args:
        df: a DataFrame
        features_ref: Features reference class.
    """
    from bionty import CellMarker, Gene, Protein

    features_ref = features_ref._entity
    df_curated = features_ref.curate(df=df, **curate_kwargs)
    if hasattr(features_ref, "_entity"):
        parsing_id = features_ref._entity._parsing_id
    else:
        parsing_id = features_ref._parsing_id

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

    features_type = (
        features_ref.entity
        if hasattr(features_ref, "entity")
        else features_ref._entity.entity
    )

    features = select(
        Features,
        id=features_hash,
        type=features_type,
    ).one_or_none()
    if features is not None:
        return features  # features already exists!
    else:
        features = Features(id=features_hash, type=features_type)
        records = get_features_records(parsing_id, features_ref, df_curated)

        if isinstance(features_ref, Gene):
            for record in records:
                features.genes.append(record)
        elif isinstance(features_ref, Protein):
            for record in records:
                features.proteins.append(record)
        elif isinstance(features_ref, CellMarker):
            for record in records:
                features.cell_markers.append(record)

    return features


def get_features(file_privates, features_ref, **curate_kwargs):
    """Updates file in place."""
    memory_rep = file_privates["_memory_rep"]
    if memory_rep is None:
        memory_rep = load_to_memory(file_privates["_local_filepath"])
    try:
        df = getattr(memory_rep, "var")  # for AnnData objects
        if callable(df):
            df = memory_rep
    except AttributeError:
        df = memory_rep
    return parse_features(df, features_ref, **curate_kwargs)
