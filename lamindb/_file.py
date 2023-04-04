from pathlib import Path
from typing import Any, Optional, Tuple, Union

import pandas as pd
from anndata import AnnData
from lamin_logger import logger
from lndb import settings as setup_settings
from lndb_storage import UPath
from lndb_storage.object import infer_suffix, size_adata, write_to_file
from lnschema_core import File, Run, Storage

from lamindb._features import get_features
from lamindb._settings import settings
from lamindb.dev.db._add import get_storage_root_and_root_str
from lamindb.dev.db._select import select
from lamindb.dev.hashing import hash_file

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
    data: Union[Path, UPath, str, pd.DataFrame, AnnData], name, format
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
        name = filepath.name
        # also see tests/test_file_hashing.py
        suffix = "".join(filepath.suffixes)
    # For now, in-memory objects are always saved to local_filepath first
    # This behavior will change in the future
    elif isinstance(data, (pd.DataFrame, AnnData)):
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
        result = select(File, hash=hash).all()
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
    memory_rep: Optional[Union[pd.DataFrame, AnnData]],
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
    data: Union[Path, UPath, str, pd.DataFrame, AnnData],
    *,
    name: Optional[str] = None,
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

    file_privates = dict(
        _local_filepath=localpath,
        _cloud_filepath=cloudpath,
        _memory_rep=memory_rep,
    )

    file_kwargs = dict(
        name=name,
        suffix=suffix,
        hash=hash,
        source_id=run.id,
        size=size,
        storage_id=storage.id,
        source=run,
    )

    return file_kwargs, file_privates


# expose to user via ln.Features
def get_features_from_data(
    data: Union[Path, UPath, str, pd.DataFrame, AnnData],
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
