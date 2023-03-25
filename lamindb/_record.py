import base64
import hashlib
from pathlib import Path
from typing import Any, List, Optional, Set, Tuple, Union

import anndata as ad
import pandas as pd
from lamin_logger import logger
from lndb import settings as setup_settings
from lndb_storage import UPath, load_to_memory
from lndb_storage.object import infer_suffix, size_adata, write_to_file
from lnschema_core import Features
from lnschema_core import File as lns_File
from lnschema_core import Run, Storage

from ._settings import settings
from .dev._core import get_name_suffix_from_filepath
from .dev.db._add import add, get_storage_root_and_root_str
from .dev.db._select import select
from .schema._table import table_meta

NO_NAME_ERROR = """
Pass a name in `ln.File(..., name=name)` when ingesting in-memory data.
"""

NO_SOURCE_ERROR = """
Error: Please link a data source using the `source` argument.
Fix: Link a data source by passing a run, e.g., via

run = lns.Run(transform=transform)
file = ln.File(..., source=run)

Or, by calling ln.context.track(), which sets a global run context.

More details: https://lamin.ai/docs/faq/ingest
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


def get_hash(local_filepath, suffix, check_hash: bool = True):
    if suffix != ".zarr":  # if not streamed
        hash = hash_file(local_filepath)
        if not check_hash:
            return hash
        result = select(lns_File, hash=hash).all()
        if len(result) > 0:
            msg = f"A file with same hash is already in the DB: {result}"
            if settings.error_on_file_hash_exists:
                hint = (
                    "💡 You can make this error a warning:\n"
                    "    ln.settings.error_on_file_hash_exists = False"
                )
                raise RuntimeError(f"{msg}\n{hint}")
            else:
                logger.warning(msg)
    else:
        hash = None
    return hash


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


def parse_features(df: pd.DataFrame, features_ref: Any, **curate_kwargs) -> None:
    """Link features to a knowledge table.

    Args:
        df: a DataFrame
        features_ref: Features reference class.
    """
    from bionty import CellMarker, Gene, Protein

    # Add and curate features against a knowledge table
    # column = None
    # if parsing_id in df.columns:
    #     column = parsing_id
    # else:
    #     logger.warning(f"{parsing_id} column not found, using index as features.")
    df_curated = features_ref.curate(df=df, **curate_kwargs)
    # TODO: fix in the next PR
    parsing_id = df_curated.index.name
    if parsing_id is None:
        if features_ref.entity == "gene":
            parsing_id = "ensembl_gene_id"
        elif features_ref.entity == "cell_marker":
            parsing_id = "name"

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
    else:
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


def get_run(run: Optional[Run]) -> Run:
    if run is None:
        from ._context import context

        run = context.run
        if run is None:
            raise ValueError(NO_SOURCE_ERROR)
    # the following ensures that queried objects (within __init__)
    # behave like queried objects, only example right now: Run
    if hasattr(run, "_ln_identity_key") and run._ln_identity_key is not None:
        run._sa_instance_state.key = run._ln_identity_key
    return run


def get_path_size_hash(
    filepath: Union[Path, UPath],
    memory_rep: Optional[Union[pd.DataFrame, ad.AnnData]],
    suffix: str,
    check_hash: bool = True,
):
    cloudpath = None
    localpath = None

    path = UPath(filepath)  # returns Path for local

    if suffix == ".zarr":
        if memory_rep is not None:
            size = size_adata(memory_rep)
        else:
            if isinstance(path, UPath):
                cloudpath = filepath
                # todo: properly calculate size
                size = 0
            else:
                localpath = filepath
                size = sum(f.stat().st_size for f in path.rglob("*") if f.is_file())
        hash = None
    else:
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
            hash = get_hash(filepath, suffix, check_hash=check_hash)

    return localpath, cloudpath, size, hash


# expose to user via ln.File
def get_file_kwargs_from_data(
    data: Union[Path, UPath, str, pd.DataFrame, ad.AnnData],
    *,
    name: Optional[str] = None,
    source: Optional[Run] = None,
    format: Optional[str] = None,
    # backward compat
    features_ref: Optional[Any] = None,
    **curate_kwargs,
):
    run = get_run(source)
    memory_rep, filepath, name, suffix = serialize(data, name, format)
    localpath, cloudpath, size, hash = get_path_size_hash(filepath, memory_rep, suffix)

    # if local_filepath is already in the configured storage location
    # skip the upload
    _, root_str = get_storage_root_and_root_str()
    storage = select(Storage, root=root_str).one()

    file_privates = dict(
        _local_filepath=localpath,
        _cloud_filepath=cloudpath,
        _memory_rep=memory_rep,
    )

    # TODO: remove later
    # backward compat
    if features_ref is not None:
        logger.warning(
            "DeprecationWarning: `features_ref` is deprecated, please use"
            " `ln.Features`!"
        )
        features = [
            get_features(file_privates, features_ref, **curate_kwargs)
        ]  # has to be list!
    else:
        features = []

    file_kwargs = dict(
        name=name,
        suffix=suffix,
        hash=hash,
        source_id=run.id,
        size=size,
        storage_id=storage.id,
        source=run,
        features=features,
    )

    return file_kwargs, file_privates


# expose to user via ln.Features
def get_features_from_data(
    data: Union[Path, UPath, str, pd.DataFrame, ad.AnnData],
    reference: Any,
    format: Optional[str] = None,
):
    memory_rep, filepath, _, suffix = serialize(data, "features", format)
    localpath, cloudpath, _, _ = get_path_size_hash(
        filepath, memory_rep, suffix, check_hash=False
    )

    file_privates = dict(
        _local_filepath=localpath,
        _cloud_filepath=cloudpath,
        _memory_rep=memory_rep,
    )
    return get_features(file_privates, reference)


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
