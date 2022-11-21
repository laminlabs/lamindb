import base64
import hashlib
from pathlib import Path
from typing import Any, List, Optional, Set, Tuple, Union

import anndata as ad
import bcoding
import lnschema_bionty as bionty
import pandas as pd
from lamin_logger import logger
from lndb_setup import settings
from lnschema_core import DObject, Features, Run, Storage
from typeguard import typechecked

from lamindb.knowledge import CellMarker, Gene, Protein

from .dev._core import get_name_suffix_from_filepath
from .dev.db._add import add
from .dev.db._select import select
from .dev.file import load_to_memory
from .dev.object import infer_suffix, size_adata, write_to_file
from .schema._table import table_meta


def serialize(
    data: Union[Path, str, pd.DataFrame, ad.AnnData], name, format
) -> Tuple[Any, Path, str, str]:
    """Serialize a data object that's provided as file or in memory."""
    memory_rep = None
    if isinstance(data, (Path, str)):
        local_filepath = Path(data)
        name, suffix = get_name_suffix_from_filepath(local_filepath)
    elif isinstance(data, (pd.DataFrame, ad.AnnData)):
        if name is None:
            raise RuntimeError("Provide name if recording in-memory data.")
        memory_rep = data
        suffix = infer_suffix(data, format)
        local_filepath = Path(f"{name}{suffix}")
        if suffix != ".zarr":
            write_to_file(data, local_filepath)
    else:
        raise NotImplementedError("Recording not yet implemented for this type.")
    return memory_rep, local_filepath, name, suffix


def get_hash(local_filepath, suffix):
    if suffix != ".zarr":  # if not streamed
        hash = hash_file(local_filepath)
        result = select(DObject, hash=hash).one_or_none()
        if result is not None:
            raise RuntimeError(
                "Based on the MD5 hash, the exact same data object is already"
                f" in the database: {result}"
            )
    else:
        hash = None
    return hash


def get_features_records(
    parsing_id: str,
    features_ref: Union[Gene, Protein, CellMarker],
    df_curated: pd.DataFrame,
) -> List[Union[Gene, Protein, CellMarker]]:
    # insert species entry if not exists
    species = select(bionty.Species, common_name=features_ref.species).one_or_none()
    if species is None:
        species = add(bionty.Species(common_name=features_ref.species))

    model = table_meta.get_model(f"bionty.{features_ref.entity}")

    # all existing feature records of the species in the db
    stmt = (
        select(model)
        .where(getattr(model, parsing_id).in_(df_curated.index))
        .where(getattr(model, "species_id") == species.id)
    )
    records = stmt.all()
    records_df = df_curated.index.intersection(stmt.df()[parsing_id])

    # new records to be appended
    new_ids = df_curated.index.difference(records_df)
    if len(new_ids) > 0:
        # mapped new_ids
        mapped = features_ref.df.loc[features_ref.df.index.intersection(new_ids)].copy()
        mapped.index.name = parsing_id
        if mapped.shape[0] > 0:
            for kwargs in mapped.reset_index().to_dict(orient="records"):
                kwargs["species_id"] = species.id
                record = model(**kwargs)
                records.append(record)
        # unmapped new_ids
        unmapped = set(new_ids).difference(mapped.index)
        if len(unmapped) > 0:
            for i in unmapped:
                record = model(**{parsing_id: i, "species_id": species.id})
                records.append(record)

    return records


def parse_features(
    df: pd.DataFrame, features_ref: Union[CellMarker, Gene, Protein]
) -> None:
    """Link features to a knowledge table.

    Args:
        df: a DataFrame
        features_ref: Features reference class.
    """
    parsing_id = features_ref._id_field

    # Add and curate features against a knowledge table
    column = None
    if parsing_id in df.columns:
        column = parsing_id
    else:
        logger.warning(f"{parsing_id} column not found, using index as features.")
    df_curated = features_ref.curate(df=df, column=column)

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

    features = select(
        Features,
        id=features_hash,
        type=features_ref.entity,
    ).one_or_none()
    if features is not None:
        return features  # features already exists!

    features = Features(id=features_hash, type=features_ref.entity)
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


def get_features(dobject, features_ref):
    """Updates dobject in place."""
    memory_rep = dobject._memory_rep
    if memory_rep is None:
        memory_rep = load_to_memory(dobject._local_filepath)
    try:
        df = getattr(memory_rep, "var")  # for AnnData objects
        if callable(df):
            df = memory_rep
    except AttributeError:
        df = memory_rep
    return parse_features(df, features_ref)


def get_run(run: Optional[Run]) -> Run:
    if run is None:
        from ._nb import _run

        run = _run
        if run is None:
            raise ValueError("Pass a Run record.")
    return run


@typechecked
def record(
    data: Union[Path, str, pd.DataFrame, ad.AnnData],
    *,
    name: Optional[str] = None,
    features_ref: Optional[Union[CellMarker, Gene, Protein]] = None,
    run: Optional[Run] = None,
    id: Optional[str] = None,
    format: Optional[str] = None,
) -> DObject:
    """Record a data object.

    Guide: :doc:`/db/guide/ingest`.

    Args:
        data: Filepath or in-memory data.
        name: Name of the data object, required if an in-memory object is passed.
        features_ref: Reference against which to link features.
        run: The data transform that links to the data source of the data object.
        id: The id of the dobject.
        format: Whether to use `h5ad` or `zarr` to store an `AnnData` object.
    """
    run = get_run(run)
    memory_rep, local_filepath, name, suffix = serialize(data, name, format)
    if suffix != ".zarr":
        size = Path(local_filepath).stat().st_size
    else:
        size = size_adata(memory_rep)
    hash = get_hash(local_filepath, suffix)
    storage = select(Storage, root=str(settings.instance.storage_root)).one()
    dobject = DObject(
        name=name,
        suffix=suffix,
        hash=hash,
        run_id=run.id,
        size=size,
        storage_id=storage.id,
        run=run,
    )
    if id is not None:  # cannot pass it into constructor due to default factory
        dobject.id = id
    dobject._local_filepath = local_filepath
    dobject._memory_rep = memory_rep
    if features_ref is not None:
        dobject.features.append(get_features(dobject, features_ref))
    return dobject


def to_b64_str(bstr: bytes):
    b64 = base64.urlsafe_b64encode(bstr).decode().strip("=")
    return b64


# a lot to read about this
# lamin-notes/2022/hashing
# redun
# bcoding
# bencode
# etc.
def hash_set(s: Set) -> str:
    return to_b64_str(hashlib.sha512(bcoding.bencode(s)).digest())[:20]


def hash_file(path: Path) -> str:
    # based on https://stackoverflow.com/questions/3431825/generating-an-md5-hash-of-a-file  # noqa
    hash_md5 = hashlib.md5()
    with open(path, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            hash_md5.update(chunk)
    return to_b64_str(hash_md5.digest())
