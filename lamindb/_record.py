import base64
import hashlib
from pathlib import Path
from typing import Any, List, Optional, Set, Tuple, Union

import anndata as ad
import lnschema_bionty as bionty
import pandas as pd
from lamin_logger import logger
from lndb import settings as setup_settings
from lndb.dev import UPath
from lnschema_core import DObject as lns_DObject
from lnschema_core import Features, Run, Storage

from lamindb.knowledge import CellMarker, Gene, Protein

from ._settings import settings
from .dev._core import get_name_suffix_from_filepath
from .dev.db._add import add, get_storage_root_and_root_str
from .dev.db._select import select
from .dev.file import load_to_memory
from .dev.object import infer_suffix, size_adata, write_to_file
from .schema._table import table_meta

NO_NAME_ERROR = """
Pass a name in `ln.DObject(..., name=name)` when ingesting in-memory data.
"""

NO_SOURCE_ERROR = """
Pass a run in `ln.DObject(..., source=run)`.
Or, if you're in a notebook, call `ln.nb.header()` at the top.
"""


def serialize(
    data: Union[Path, UPath, str, pd.DataFrame, ad.AnnData], name, format
) -> Tuple[Any, Union[Path, UPath], str, str]:
    """Serialize a data object that's provided as file or in memory."""
    # Convert str to either Path or CloudPath
    if isinstance(data, (str, Path, UPath)):
        filepath = UPath(data)  # returns Path for local
        if (
            isinstance(filepath, UPath)
            and setup_settings.instance.storage.root not in filepath.parents  # noqa
        ):
            raise ValueError(
                "Can only track objects in configured cloud storage locations."
                " Please call `lndb.set_storage('< bucket_name >')`."
            )
        memory_rep = None
        name, suffix = get_name_suffix_from_filepath(filepath)
    # For now, in-memory objects are always saved to local_filepath first
    # This behavior will change in the future
    elif isinstance(data, (pd.DataFrame, ad.AnnData)):
        if name is None:
            raise ValueError(NO_NAME_ERROR)
        memory_rep = data
        suffix = infer_suffix(data, format)
        # this is always local
        filepath = Path(f"{name}{suffix}")
        if suffix != ".zarr":
            write_to_file(data, filepath)
    else:
        raise NotImplementedError("Recording not yet implemented for this type.")
    return memory_rep, filepath, name, suffix


def get_hash(local_filepath, suffix):
    if suffix != ".zarr":  # if not streamed
        hash = hash_file(local_filepath)
        result = select(lns_DObject, hash=hash).all()
        if len(result) > 0:
            msg = f"A dobject with same hash is already in the DB: {result}"
            if settings.error_on_dobject_hash_exists:
                hint = (
                    "💡 You can make this error a warning:\n"
                    "    ln.settings.error_on_dobject_hash_exists = False"
                )
                raise RuntimeError(f"{msg}\n{hint}")
            else:
                logger.warning(msg)
    else:
        hash = None
    return hash


def get_features_records(
    parsing_id: str,
    features_ref: Union[Gene, Protein, CellMarker],
    df_curated: pd.DataFrame,
) -> List[Union[Gene, Protein, CellMarker]]:
    # insert species entry if not exists
    species = add(bionty.Species, name=features_ref.species)

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
        mapped = features_ref.df.loc[features_ref.df.index.intersection(new_ids)].copy()
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


def get_features(dobject_privates, features_ref):
    """Updates dobject in place."""
    memory_rep = dobject_privates["_memory_rep"]
    if memory_rep is None:
        memory_rep = load_to_memory(dobject_privates["_local_filepath"])
    try:
        df = getattr(memory_rep, "var")  # for AnnData objects
        if callable(df):
            df = memory_rep
    except AttributeError:
        df = memory_rep
    return parse_features(df, features_ref)


def get_run(run: Optional[Run]) -> Run:
    if run is None:
        from . import nb

        run = nb.run
        if run is None:
            raise ValueError(NO_SOURCE_ERROR)
    return run


def get_path_size_hash(
    filepath: Union[Path, UPath],
    memory_rep: Optional[Union[pd.DataFrame, ad.AnnData]],
    suffix: str,
):
    cloudpath = None
    localpath = None
    # both cloudpath and localpath are None for zarr
    if suffix != ".zarr":
        path = UPath(filepath)  # returns Path for local
        if isinstance(path, UPath):
            try:
                size = path.stat()["size"]
            # here trying to fix access issue with new s3 buckets
            except Exception as e:
                if path._url.scheme == "s3":
                    path = UPath(filepath, cache_regions=True)
                    size = path.stat()["size"]
                else:
                    raise e
            cloudpath = filepath
            hash = None
        else:
            size = path.stat().st_size
            localpath = filepath
            hash = get_hash(filepath, suffix)
    else:
        size = size_adata(memory_rep)
        hash = None

    return localpath, cloudpath, size, hash


def get_dobject_kwargs_from_data(
    data: Union[Path, UPath, str, pd.DataFrame, ad.AnnData],
    *,
    name: Optional[str] = None,
    features_ref: Optional[Union[CellMarker, Gene, Protein]] = None,
    source: Optional[Run] = None,
    format: Optional[str] = None,
):
    run = get_run(source)
    memory_rep, filepath, name, suffix = serialize(data, name, format)
    localpath, cloudpath, size, hash = get_path_size_hash(filepath, memory_rep, suffix)

    # if local_filepath is already in the configured storage location
    # skip the upload
    _, root_str = get_storage_root_and_root_str()
    storage = select(Storage, root=root_str).one()

    dobject_privates = dict(
        _local_filepath=localpath,
        _cloud_filepath=cloudpath,
        _memory_rep=memory_rep,
    )

    if features_ref is not None:
        features = [get_features(dobject_privates, features_ref)]  # has to be list!
    else:
        features = []
    dobject_kwargs = dict(
        name=name,
        suffix=suffix,
        hash=hash,
        source_id=run.id,
        size=size,
        storage_id=storage.id,
        source=run,
        features=features,
    )
    return dobject_kwargs, dobject_privates


def to_b64_str(bstr: bytes):
    b64 = base64.urlsafe_b64encode(bstr).decode().strip("=")
    return b64


# a lot to read about this: lamin-notes/2022/hashing
def hash_set(s: Set[str]) -> str:
    bstr = ":".join(sorted(s)).encode("utf-8")
    # as we're truncating at 20 b64, we choose md5 over sha512
    return to_b64_str(hashlib.md5(bstr).digest())[:20]


def hash_file(path: Path) -> str:
    # based on https://stackoverflow.com/questions/3431825/generating-an-md5-hash-of-a-file  # noqa
    hash_md5 = hashlib.md5()
    with open(path, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            hash_md5.update(chunk)
    return to_b64_str(hash_md5.digest())
