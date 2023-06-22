from pathlib import Path, PurePath, PurePosixPath
from typing import Any, Optional, Tuple, Union

import lamindb_setup
import pandas as pd
from anndata import AnnData
from appdirs import AppDirs
from lamin_logger import logger
from lnschema_core import File, Run, ids
from lnschema_core.types import DataLike, PathLike

from lamindb._file_access import auto_storage_key_from_file
from lamindb.dev._settings import settings
from lamindb.dev.hashing import hash_file
from lamindb.dev.storage import UPath, infer_suffix, size_adata, write_to_file

from ._file_access import AUTO_KEY_PREFIX

DIRS = AppDirs("lamindb", "laminlabs")

NO_NAME_ERROR = """\
Pass a name or key in `ln.File(...)`.
"""


def serialize(
    data: Union[Path, UPath, str, pd.DataFrame, AnnData],
    name: Optional[str],
    format,
    key: Optional[str],
) -> Tuple[Any, Union[Path, UPath], str]:
    """Serialize a data object that's provided as file or in memory."""
    # Convert str to either Path or UPath
    if isinstance(data, (str, Path, UPath)):
        filepath = UPath(data)  # returns Path for local
        try:  # check if file exists
            if not filepath.exists():
                raise FileNotFoundError(filepath)
        except PermissionError:  # we will setup permissions later
            pass
        if isinstance(filepath, UPath):
            new_storage = list(filepath.parents)[-1]
            if not get_check_path_in_storage(filepath):
                raise ValueError(
                    "Currently do not support moving cloud data across buckets."
                    " Set default storage to point to your cloud bucket:\n"
                    f" `ln.setup.set.storage({new_storage})` or `lamin set --storage"
                    f" {new_storage}`"
                )
            root = lamindb_setup.settings.storage.root
            if isinstance(root, UPath):
                filepath = UPath(
                    filepath, **root._kwargs
                )  # inherit fsspec kwargs from root
        memory_rep = None
        # also see tests/test_file_hashing.py
        suffix = "".join(filepath.suffixes)
    # For now, in-memory objects are always saved to local_filepath first
    # This behavior will change in the future
    elif isinstance(data, (pd.DataFrame, AnnData)):
        memory_rep = data
        suffix = infer_suffix(data, format)
        # the following filepath is always local
        if name is None and key is None:
            raise ValueError(NO_NAME_ERROR)
        elif name is not None:
            cache_name = name
        else:
            cache_name = Path(key).name
        if lamindb_setup.settings.storage.cache_dir is not None:
            filepath = lamindb_setup.settings.storage.cache_dir / cache_name
        else:
            # this should likely be added to lamindb_setup.settings.storage
            cache_dir = Path(DIRS.user_cache_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)
            filepath = cache_dir / cache_name
        if filepath.suffixes == []:
            filepath = filepath.with_suffix(suffix)
        if suffix != ".zarr":
            write_to_file(data, filepath)
    else:
        raise NotImplementedError("Recording not yet implemented for this type.")
    return memory_rep, filepath, suffix


def get_hash(
    local_filepath, suffix, check_hash: bool = True
) -> Optional[Union[str, File]]:
    if suffix in {".zarr", ".zrad"}:
        return None
    hash = hash_file(local_filepath)
    if not check_hash:
        return hash
    result = File.select(hash=hash).list()
    if len(result) > 0:
        if settings.upon_file_create_if_hash_exists == "error":
            msg = f"A file with same hash exists: {result[0]}"
            hint = (
                "💡 You can make this error a warning:\n"
                "    ln.settings.upon_file_create_if_hash_exists"
            )
            raise RuntimeError(f"{msg}\n{hint}")
        elif settings.upon_file_create_if_hash_exists == "warn_create_new":
            logger.warning(
                "Creating new File object despite existing file with same hash:"
                f" {result[0]}"
            )
            return hash
        else:
            logger.warning(f"Returning existing file with same hash: {result[0]}")
            return result[0]
    else:
        return hash


def get_run(run: Optional[Run]) -> Optional[Run]:
    if run is None:
        from ._context import context

        run = context.run
        if run is None:
            logger.hint(
                "no run & transform get linked, consider passing a `run` or calling"
                " ln.track()"
            )
    return run


def get_path_size_hash(
    filepath: Union[Path, UPath],
    memory_rep: Optional[Union[pd.DataFrame, AnnData]],
    suffix: str,
    check_hash: bool = True,
):
    cloudpath = None
    localpath = None

    if suffix == ".zarr":
        if memory_rep is not None:
            size = size_adata(memory_rep)
        else:
            if isinstance(filepath, UPath):
                cloudpath = filepath
                # todo: properly calculate size
                size = 0
            else:
                localpath = filepath
                size = sum(
                    f.stat().st_size for f in filepath.rglob("*") if f.is_file()  # type: ignore # noqa
                )
        hash = None
    else:
        if isinstance(filepath, UPath):
            try:
                size = filepath.stat()["size"]
            # here trying to fix access issue with new s3 buckets
            except Exception as e:
                if filepath._url.scheme == "s3":
                    filepath = UPath(filepath, cache_regions=True)
                    size = filepath.stat()["size"]
                else:
                    raise e
            cloudpath = filepath
            hash = None
        else:
            size = filepath.stat().st_size  # type: ignore
            localpath = filepath
            hash = get_hash(filepath, suffix, check_hash=check_hash)

    return localpath, cloudpath, size, hash


def get_check_path_in_storage(
    filepath: Union[Path, UPath], *, root: Optional[Union[Path, UPath]] = None
) -> bool:
    assert isinstance(filepath, Path)
    if root is None:
        root = lamindb_setup.settings.storage.root
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
        root = lamindb_setup.settings.storage.root
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
    memory_rep, filepath, suffix = serialize(data, name, format, key)
    # the following will return a localpath that is not None if filepath is local
    # it will return a cloudpath that is not None if filepath is on the cloud
    local_filepath, cloud_filepath, size, hash = get_path_size_hash(
        filepath, memory_rep, suffix
    )
    if isinstance(hash, File):
        return hash, None
    check_path_in_storage = get_check_path_in_storage(filepath)
    # if we pass a file, no storage key, and path is already in storage,
    # then use the existing relative path within the storage location
    # as storage key
    if memory_rep is None and key is None and check_path_in_storage:
        key = get_relative_path_to_root(path=filepath).as_posix()

    if name is None and key is None:
        raise ValueError(NO_NAME_ERROR)

    if key is not None and key.startswith(AUTO_KEY_PREFIX):
        raise ValueError(f"Key cannot start with {AUTO_KEY_PREFIX}")

    kwargs = dict(
        name=name,
        suffix=suffix,
        hash=hash,
        key=key,
        size=size,
        storage_id=lamindb_setup.settings.storage.id,
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


def log_storage_hint(
    *, check_path_in_storage: bool, key: str, id: str, suffix: str
) -> None:
    hint = ""
    if check_path_in_storage:
        hint += "file in storage ✓"
    else:
        hint += "file will be copied to storage upon `save()`"
    if key is None:
        hint += f" using storage key = {id}{suffix}"
    else:
        hint += f" using storage key = {key}"
    logger.hint(hint)


def init_file(file: File, *args, **kwargs):
    # Below checks for the Django-internal call in from_db()
    # it'd be better if we could avoid this, but not being able to create a File
    # from data with the default constructor renders the central class of the API
    # essentially useless
    # The danger below is not that a user might pass as many args (12 of it), but rather
    # that at some point the Django API might change; on the other hand, this
    # condition of for calling the constructor based on kwargs should always
    # stay robust
    if len(args) == len(file._meta.concrete_fields):
        super(File, file).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 1:
        raise ValueError("Only one non-keyword arg allowed: data")
    data: Union[PathLike, DataLike] = kwargs["data"] if len(args) == 0 else args[0]
    key: Optional[str] = kwargs["key"] if "key" in kwargs else None
    name: Optional[str] = kwargs["name"] if "name" in kwargs else None
    run: Optional[Run] = kwargs["run"] if "run" in kwargs else None
    format = kwargs["format"] if "format" in kwargs else None

    kwargs, privates = get_file_kwargs_from_data(
        data=data,
        name=name,
        key=key,
        run=run,
        format=format,
    )
    # an object with the same hash already exists
    if isinstance(kwargs, File):
        # this is the way Django instantiates from the DB internally
        # https://github.com/django/django/blob/549d6ffeb6d626b023acc40c3bb2093b4b25b3d6/django/db/models/base.py#LL488C1-L491C51
        new_args = [
            getattr(kwargs, field.attname) for field in file._meta.concrete_fields
        ]
        super(File, file).__init__(*new_args)
        return None

    kwargs["id"] = ids.base62_20()
    log_storage_hint(
        check_path_in_storage=privates["check_path_in_storage"],
        key=kwargs["key"],
        id=kwargs["id"],
        suffix=kwargs["suffix"],
    )

    # transform cannot be directly passed, just via run
    # it's directly stored in the file table to avoid another join
    # mediate by the run table
    if kwargs["run"] is not None:
        if kwargs["run"].transform_id is not None:
            kwargs["transform_id"] = kwargs["run"].transform_id
        else:
            # accessing the relationship should always be possible if
            # the above if clause was false as then, we should have a fresh
            # Transform object that is not queried from the DB
            assert kwargs["run"].transform is not None
            kwargs["transform"] = kwargs["run"].transform

    if data is not None:
        file._local_filepath = privates["local_filepath"]
        file._cloud_filepath = privates["cloud_filepath"]
        file._memory_rep = privates["memory_rep"]
        file._to_store = not privates["check_path_in_storage"]

    super(File, file).__init__(**kwargs)


def from_dir(
    path: Union[Path, UPath, str],
    *,
    run: Optional[Run] = None,
):
    """Create file records from a directory."""
    folderpath = UPath(path)
    check_path_in_storage = get_check_path_in_storage(folderpath)

    if check_path_in_storage:
        folder_key = get_relative_path_to_root(path=folderpath).as_posix()
    else:
        raise RuntimeError(
            "Currently, only directories in default storage can be registered!\n"
            "You can either move your folder into the current default storage"
            "or add a new default storage through `ln.settings.storage`"
        )
    # always sanitize by stripping a trailing slash
    folder_key = folder_key.rstrip("/")
    logger.hint(f"using storage prefix = {folder_key}/")

    # TODO: UPath doesn't list the first level files and dirs with "*"
    pattern = "" if isinstance(folderpath, UPath) else "*"

    # silence fine-grained logging
    verbosity = settings.verbosity
    settings.verbosity = 1  # just warnings
    files = []
    for filepath in folderpath.rglob(pattern):
        if filepath.is_file():
            relative_path = get_relative_path_to_directory(filepath, folderpath)
            file_key = folder_key + "/" + relative_path.as_posix()
            files.append(File(filepath, run=run, key=file_key))
    settings.verbosity = verbosity
    logger.info(f"→ {len(files)} files")
    return files


def replace_file(
    file: File,
    data: Union[PathLike, DataLike] = None,
    run: Optional[Run] = None,
    format: Optional[str] = None,
):
    kwargs, privates = get_file_kwargs_from_data(
        data=data,
        name=file.name,
        key=file.key,
        run=run,
        format=format,
    )
    if file.key is not None:
        key_path = PurePosixPath(file.key)
        if isinstance(data, (Path, str)) and kwargs["name"] is not None:
            new_name = kwargs["name"]  # use the name from the data filepath
        else:
            # do not change the key stem to file.name
            new_name = key_path.stem  # use the stem of the key for in-memory data
        if PurePosixPath(new_name).suffixes == []:
            new_name = f"{new_name}{kwargs['suffix']}"
        if key_path.name != new_name:
            file._clear_storagekey = file.key
            file.key = str(key_path.with_name(new_name))
            logger.warning(
                f"Replacing the file will replace key '{key_path}' with '{file.key}'"
                f" and delete '{key_path}' upon `save()`"
            )
    else:
        file.key = kwargs["key"]
        old_storage = auto_storage_key_from_file(file)
        new_storage = (
            file.key if file.key is not None else f"{file.id}{kwargs['suffix']}"
        )
        if old_storage != new_storage:
            file._clear_storagekey = old_storage

    file.suffix = kwargs["suffix"]
    file.size = kwargs["size"]
    file.hash = kwargs["hash"]
    file.run = kwargs["run"]
    file._local_filepath = privates["local_filepath"]
    file._cloud_filepath = privates["cloud_filepath"]
    file._memory_rep = privates["memory_rep"]
    file._to_store = not privates[
        "check_path_in_storage"
    ]  # no need to upload if new file is already in storage
