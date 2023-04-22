from pathlib import Path, PurePath
from typing import Any, Optional, Tuple, Union

import lndb
import pandas as pd
from anndata import AnnData
from appdirs import AppDirs
from lamin_logger import logger
from lndb_storage import UPath
from lndb_storage.object import infer_suffix, size_adata, write_to_file
from lnschema_core import File, Run

from lamindb._features import get_features
from lamindb._settings import settings
from lamindb.dev.db._select import select
from lamindb.dev.hashing import hash_file

DIRS = AppDirs("lamindb", "laminlabs")

NO_NAME_ERROR = """
Pass a name or key in `ln.File(...)` when ingesting in-memory data.
"""

NO_SOURCE_HINT = """
Consider linking a run using the `run` argument, e.g., by calling ln.track().
"""


def serialize(
    data: Union[Path, UPath, str, pd.DataFrame, AnnData],
    name,
    format,
    key: Optional[str],
) -> Tuple[Any, Union[Path, UPath], str, str]:
    """Serialize a data object that's provided as file or in memory."""
    # Convert str to either Path or CloudPath
    if isinstance(data, (str, Path, UPath)):
        filepath = UPath(data)  # returns Path for local
        if (
            isinstance(filepath, UPath)
            and lndb.settings.instance.storage.root not in filepath.parents  # noqa
        ):
            raise ValueError(
                "Can only track objects in configured cloud storage locations."
                " Please call `lndb.set_storage('< bucket_name >')`."
            )
        memory_rep = None
        if key is None:
            name = filepath.name
        else:
            name = PurePath(key).name
        # also see tests/test_file_hashing.py
        suffix = "".join(filepath.suffixes)
    # For now, in-memory objects are always saved to local_filepath first
    # This behavior will change in the future
    elif isinstance(data, (pd.DataFrame, AnnData)):
        if name is None and key is None:
            raise ValueError(NO_NAME_ERROR)
        if name is None:
            name = PurePath(key).name
        memory_rep = data
        suffix = infer_suffix(data, format)
        # the following filepath is always local
        if lndb.settings.storage.cache_dir is not None:
            filepath = lndb.settings.storage.cache_dir / name
        else:
            # this should likely be added to lndb.settings.storage
            cache_dir = Path(DIRS.user_cache_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)
            filepath = cache_dir / name
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


def get_run(run: Optional[Run]) -> Optional[Run]:
    if run is None:
        from ._context import context

        run = context.run
        if run is None:
            logger.warning("No run & transform get linked to this file.")
            logger.hint(NO_SOURCE_HINT)
    # the following ensures that queried objects (within __init__)
    # behave like queried objects, only example right now: Run
    if hasattr(run, "_ln_identity_key") and run._ln_identity_key is not None:  # type: ignore  # noqa
        run._sa_instance_state.key = run._ln_identity_key  # type: ignore
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


def get_check_path_in_storage(
    filepath: Union[Path, UPath], *, root: Optional[Union[Path, UPath]] = None
) -> bool:
    assert isinstance(filepath, Path)
    if root is None:
        root = lndb.settings.storage.root
    # the following comparisons can fail if types aren't comparable
    if isinstance(filepath, UPath) and isinstance(root, UPath):
        # the following tests equivalency of two UPath objects
        # via string representations; otherwise
        # S3Path('s3://lndb-storage/') and S3Path('s3://lamindb-ci/')
        # test as equivalent
        return list(filepath.parents)[-1].as_posix() == root.as_posix()
    elif not isinstance(filepath, UPath) and not isinstance(root, UPath):
        return root in filepath.resolve().parents
    else:
        return False


def get_relative_path_to_directory(
    path: Union[PurePath, Path, UPath], directory: Union[PurePath, Path, UPath]
) -> Union[PurePath, Path]:
    if isinstance(directory, UPath):
        # UPath.relative_to() is not behaving as it should (2023-04-07)
        relpath = PurePath(path.as_posix().replace(directory.as_posix(), ""))
    elif isinstance(directory, Path):
        relpath = path.resolve().relative_to(directory.resolve())  # type: ignore
    elif isinstance(directory, PurePath):
        relpath = path.relative_to(directory)
    else:
        raise TypeError("directory not of type Path or UPath")
    return relpath


def get_relative_path_to_root(
    path: Union[Path, UPath], *, root: Optional[Union[Path, UPath]] = None
) -> Union[PurePath, Path]:
    """Relative path to the storage root path."""
    if root is None:
        root = lndb.settings.storage.root
    return get_relative_path_to_directory(path, root)


# expose to user via ln.File
def get_file_kwargs_from_data(
    data: Union[Path, UPath, str, pd.DataFrame, AnnData],
    *,
    name: Optional[str] = None,
    key: Optional[str] = None,
    run: Optional[Run] = None,
    format: Optional[str] = None,
):
    run = get_run(run)
    memory_rep, filepath, safe_name, suffix = serialize(data, name, format, key)
    # the following will return a localpath that is not None if filepath is local
    # it will return a cloudpath that is not None if filepath is on the cloud
    local_filepath, cloud_filepath, size, hash = get_path_size_hash(
        filepath, memory_rep, suffix
    )
    check_path_in_storage = get_check_path_in_storage(filepath)
    if cloud_filepath is not None and not check_path_in_storage:
        new_storage = list(cloud_filepath.parents)[-1]
        raise ValueError(
            "Currently do not support moving cloud data across buckets. Configure"
            " storage to point to your cloud bucket:\n"
            f" `ln.setup.set.storage({new_storage})` or `lamin set --storage"
            f" {new_storage}`"
        )

    # if we pass a file, no storage key, and path is already in storage,
    # then use the existing relative path within the storage location
    # as storage key
    if memory_rep is None and key is None and check_path_in_storage:
        key = get_relative_path_to_root(path=filepath).as_posix()

    kwargs = dict(
        name=safe_name,
        suffix=suffix,
        hash=hash,
        key=key,
        size=size,
        storage_id=lndb.settings.storage.id,
        # passing both the id and the object
        # to make them both available immediately
        # after object creation
        run_id=run.id if run is not None else None,
        run=run,
    )
    privates = dict(
        local_filepath=local_filepath,
        cloud_filepath=cloud_filepath,
        memory_rep=memory_rep,
        check_path_in_storage=check_path_in_storage,
    )

    return kwargs, privates


# expose to user via ln.Features
def get_features_from_data(
    data: Union[Path, UPath, str, pd.DataFrame, AnnData],
    reference: Any,
    format: Optional[str] = None,
    **curate_kwargs,
):
    memory_rep, filepath, _, suffix = serialize(data, "features", format, key=None)
    localpath, cloudpath, _, _ = get_path_size_hash(
        filepath, memory_rep, suffix, check_hash=False
    )

    file_privates = dict(
        _local_filepath=localpath,
        _cloud_filepath=cloudpath,
        _memory_rep=memory_rep,
    )
    return get_features(file_privates, reference, **curate_kwargs)
