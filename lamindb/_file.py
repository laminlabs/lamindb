from itertools import islice
from pathlib import Path, PurePath, PurePosixPath
from typing import Any, List, Optional, Set, Tuple, Union

import anndata as ad
import lamindb_setup
import pandas as pd
from anndata import AnnData
from appdirs import AppDirs
from lamin_utils import colors, logger
from lamindb_setup import settings as setup_settings
from lamindb_setup._init_instance import register_storage
from lamindb_setup.dev import StorageSettings
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Feature, FeatureSet, File, Run, Storage, ids
from lnschema_core.types import AnnDataLike, DataLike, FieldAttr, PathLike

from lamindb.dev import run_context
from lamindb.dev._settings import settings
from lamindb.dev.hashing import b16_to_b64, hash_file
from lamindb.dev.storage import (
    UPath,
    delete_storage,
    infer_suffix,
    load_to_memory,
    size_adata,
    write_to_file,
)
from lamindb.dev.storage._backed_access import AnnDataAccessor, BackedAccessor
from lamindb.dev.storage.file import (
    ProgressCallback,
    _str_to_path,
    auto_storage_key_from_file,
    extract_suffix_from_path,
    filepath_from_file,
)
from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._feature import convert_numpy_dtype_to_lamin_feature_type
from ._parents import view_lineage
from .dev.storage.file import AUTO_KEY_PREFIX

DIRS = AppDirs("lamindb", "laminlabs")


def process_pathlike(
    filepath: Union[Path, UPath], skip_existence_check: bool = False
) -> Tuple[Storage, bool]:
    if not skip_existence_check:
        try:  # check if file exists
            if not filepath.exists():
                raise FileNotFoundError(filepath)
        except PermissionError:
            pass
    if not isinstance(filepath, UPath):
        filepath = filepath.resolve()
    # check whether the path is in default storage
    default_storage = lamindb_setup.settings.storage.record
    if check_path_is_child_of_root(filepath, default_storage.root_as_path()):
        use_existing_storage_key = True
        return default_storage, use_existing_storage_key
    else:
        # check whether the path is part of one of the existing
        # already-registered storage locations
        result = check_path_in_existing_storage(filepath)
        if isinstance(result, Storage):
            use_existing_storage_key = True
            return result, use_existing_storage_key
        else:
            # if the path is in the cloud, we have a good candidate
            # for the storage root: the bucket
            if isinstance(filepath, UPath):
                # for a UPath, new_root is always the bucket name
                new_root = list(filepath.parents)[-1]
                new_root_str = new_root.as_posix()
                logger.warning(
                    f"creating new storage location for root: {new_root_str}"
                )
                storage_settings = StorageSettings(new_root_str)
                register_storage(storage_settings)
                use_existing_storage_key = True
                return storage_settings.record, use_existing_storage_key
            # if the filepath is local
            else:
                use_existing_storage_key = False
                # if the default storage is local we'll throw an error if the user
                # doesn't provide a key
                if not lamindb_setup.settings.storage.is_cloud:
                    return default_storage, use_existing_storage_key
                # if the default storage is in the cloud (the file is going to
                # be uploaded upon saving it), we treat the filepath as a cache
                else:
                    return default_storage, use_existing_storage_key


def process_data(
    provisional_id: str,
    data: Union[PathLike, DataLike],
    format,
    skip_existence_check: bool = False,
) -> Tuple[Any, Union[Path, UPath], str, Storage, bool]:
    """Serialize a data object that's provided as file or in memory."""
    # if not overwritten, data gets stored in default storage
    if isinstance(data, (str, Path, UPath)):  # PathLike, spelled out
        filepath = (
            data if isinstance(data, (Path, UPath)) else _str_to_path(data)
        )  # returns Path for local
        storage, use_existing_storage_key = process_pathlike(
            filepath, skip_existence_check=skip_existence_check
        )
        suffix = extract_suffix_from_path(filepath)
        memory_rep = None
    elif isinstance(data, (pd.DataFrame, AnnData)):  # DataLike, spelled out
        storage = lamindb_setup.settings.storage.record
        memory_rep = data
        suffix = infer_suffix(data, format)
        cache_name = f"{provisional_id}{suffix}"
        if lamindb_setup.settings.storage.cache_dir is not None:
            filepath = lamindb_setup.settings.storage.cache_dir / cache_name
        else:
            # this should likely be added to lamindb_setup.settings.storage
            cache_dir = Path(DIRS.user_cache_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)
            filepath = cache_dir / cache_name
        # Alex: I don't understand the line below
        if filepath.suffixes == []:
            filepath = filepath.with_suffix(suffix)
        if suffix != ".zarr":
            write_to_file(data, filepath)
        use_existing_storage_key = False
    else:
        raise NotImplementedError(
            f"Do not know how to create a file object from {data}, pass a filepath"
            " instead!"
        )
    return memory_rep, filepath, suffix, storage, use_existing_storage_key


def get_hash(
    filepath,
    suffix,
    filepath_stat=None,
    check_hash: bool = True,
) -> Union[Tuple[Optional[str], Optional[str]], File]:
    if suffix in {".zarr", ".zrad"}:
        return None
    if isinstance(filepath, UPath):
        stat = filepath_stat
        if stat is not None and "ETag" in stat:
            # small files
            if "-" not in stat["ETag"]:
                # only store hash for non-multipart uploads
                # we can't rapidly validate multi-part uploaded files client-side
                # we can add more logic later down-the-road
                hash = b16_to_b64(stat["ETag"])
                hash_type = "md5"
            else:
                stripped_etag, suffix = stat["ETag"].split("-")
                suffix = suffix.strip('"')
                hash = f"{b16_to_b64(stripped_etag)}-{suffix}"
                hash_type = "md5-n"  # this is the S3 chunk-hashing strategy
        else:
            logger.warning(f"did not add hash for {filepath}")
            return None, None
    else:
        hash, hash_type = hash_file(filepath)
    if not check_hash:
        return hash, hash_type
    result = File.filter(hash=hash).list()
    if len(result) > 0:
        if settings.upon_file_create_if_hash_exists == "error":
            msg = f"file with same hash exists: {result[0]}"
            hint = (
                "💡 you can make this error a warning:\n"
                "    ln.settings.upon_file_create_if_hash_exists"
            )
            raise RuntimeError(f"{msg}\n{hint}")
        elif settings.upon_file_create_if_hash_exists == "warn_create_new":
            logger.warning(
                "creating new File object despite existing file with same hash:"
                f" {result[0]}"
            )
            return hash, hash_type
        else:
            logger.warning(f"returning existing file with same hash: {result[0]}")
            return result[0]
    else:
        return hash, hash_type


def get_run(run: Optional[Run]) -> Optional[Run]:
    if run is None:
        run = run_context.run
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
    hash_and_type: Tuple[Optional[str], Optional[str]]

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
        hash_and_type = None, None
    else:
        # to accelerate ingesting high numbers of files
        if settings.upon_file_create_skip_size_hash:
            size = None
            hash_and_type = None, None
        else:
            filepath_stat = filepath.stat()
            if isinstance(filepath, UPath):
                try:
                    size = filepath_stat["size"]  # type: ignore
                # here trying to fix access issue with new s3 buckets
                except Exception as e:
                    if filepath._url.scheme == "s3":
                        filepath = UPath(filepath, cache_regions=True)
                        size = filepath_stat["size"]  # type: ignore
                    else:
                        raise e
                cloudpath = filepath
                hash_and_type = None, None
            else:
                size = filepath_stat.st_size  # type: ignore
                localpath = filepath
            hash_and_type = get_hash(
                filepath, suffix, filepath_stat=filepath_stat, check_hash=check_hash
            )
    return localpath, cloudpath, size, hash_and_type


def check_path_in_existing_storage(
    filepath: Union[Path, UPath]
) -> Union[Storage, bool]:
    for storage in Storage.filter().all():
        # if path is part of storage, return it
        if check_path_is_child_of_root(filepath, root=_str_to_path(storage.root)):
            return storage
    return False


def check_path_is_child_of_root(
    filepath: Union[Path, UPath], root: Optional[Union[Path, UPath]] = None
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
        return root.resolve() in filepath.resolve().parents
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
        raise TypeError("Directory not of type Path or UPath")
    return relpath


def get_file_kwargs_from_data(
    *,
    data: Union[Path, UPath, str, pd.DataFrame, AnnData],
    key: Optional[str],
    run: Optional[Run],
    format: Optional[str],
    provisional_id: str,
    skip_check_exists: bool = False,
):
    run = get_run(run)
    memory_rep, filepath, suffix, storage, use_existing_storage_key = process_data(
        provisional_id, data, format, skip_check_exists
    )
    # the following will return a localpath that is not None if filepath is local
    # it will return a cloudpath that is not None if filepath is on the cloud
    local_filepath, cloud_filepath, size, hash_and_type = get_path_size_hash(
        filepath,
        memory_rep,
        suffix,
    )
    if isinstance(hash_and_type, File):
        return hash_and_type, None
    else:
        hash, hash_type = hash_and_type

    check_path_in_storage = False
    if key is None and use_existing_storage_key:
        key = get_relative_path_to_directory(
            path=filepath, directory=storage.root_as_path()
        ).as_posix()
        check_path_in_storage = True
    else:
        storage = lamindb_setup.settings.storage.record

    if key is not None and key.startswith(AUTO_KEY_PREFIX):
        raise ValueError(f"Key cannot start with {AUTO_KEY_PREFIX}")

    log_storage_hint(
        check_path_in_storage=check_path_in_storage,
        storage=storage,
        key=key,
        id=provisional_id,
        suffix=suffix,
    )

    kwargs = dict(
        suffix=suffix,
        hash=hash,
        hash_type=hash_type,
        key=key,
        size=size,
        storage_id=storage.id,
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
    *,
    check_path_in_storage: bool,
    storage: Optional[Storage],
    key: Optional[str],
    id: str,
    suffix: str,
) -> None:
    hint = ""
    if check_path_in_storage:
        hint += f"file in storage '{storage.root}'"  # type: ignore
    else:
        hint += "file will be copied to default storage upon `save()`"
    if key is None:
        hint += f" with key '{id}{suffix}'"
    else:
        hint += f" with key '{key}'"
    logger.hint(hint)


def data_is_anndata(data: DataLike):
    if isinstance(data, AnnData):
        return True
    if isinstance(data, (str, Path, UPath)):
        return Path(data).suffix in {".h5ad", ".zrad"}
    return False  # pragma: no cover


def data_is_mudata(data: DataLike):  # pragma: no cover
    try:
        from mudata import MuData
    except ModuleNotFoundError:
        return False

    if isinstance(data, MuData):
        return True
    if isinstance(data, (str, Path, UPath)):
        return Path(data).suffix in {".h5mu"}
    return False


# uses `initial_version_id` to extract a stem_id that's part of id
# this entire piece of logic might be removed in the future if it doesn't turn out
# to be robustly maintainable
def init_id(
    *,
    provisional_id: Optional[None] = None,
    initial_version_id: Optional[str] = None,
    version: Optional[str] = None,
) -> str:
    if version is not None:
        if not isinstance(version, str):
            raise ValueError(
                "`version` parameter must be `None` or `str`, e.g., '0', '1', etc."
            )
    if initial_version_id is not None:
        stem_id = initial_version_id[:18]
    else:
        stem_id = None
    # first consider an unversioned record
    if version is None and stem_id is None:
        provisional_id = ids.base62_20()
        return provisional_id  # type: ignore
    # now consider a versioned record
    id_ext = ids.base62(2)
    if provisional_id is None and stem_id is None:
        stem_id = ids.base62_18()
        provisional_id = stem_id + id_ext
    elif stem_id is not None:
        assert isinstance(stem_id, str) and len(stem_id) == 18
        provisional_id = stem_id + id_ext
    elif provisional_id is not None:
        assert isinstance(provisional_id, str) and len(provisional_id) == 20
        stem_id = provisional_id[:18]
    return provisional_id  # type: ignore


def set_version(version: Optional[str] = None, previous_version: Optional[str] = None):
    """(Auto-) set version.

    If `version` is `None`, returns the stored version.
    Otherwise sets the version to the passed version.

    Args:
        version: Version string.
        stored_version: Mock stored version for testing purposes.
    """
    if version is None and previous_version is not None:
        try:
            version = str(int(previous_version) + 1)  # increment version by 1
        except ValueError:
            raise ValueError(
                "Cannot auto-increment non-integer castable version, please provide"
                " manually"
            )
    return version


def get_ids_from_old_version_of_file(
    is_new_version_of: File,
    version: Optional[str],
) -> Tuple[str, str, str]:
    """{}"""
    msg = ""
    if is_new_version_of.version is None:
        previous_version = "1"
        msg = "setting version of old file to '1'"
    else:
        previous_version = is_new_version_of.version
    version = set_version(version, previous_version)
    if is_new_version_of.initial_version_id is None:
        initial_version_id = is_new_version_of.id
    else:
        initial_version_id = is_new_version_of.initial_version_id
    new_file_id = init_id(
        provisional_id=is_new_version_of.id,
        initial_version_id=initial_version_id,
        version=version,
    )
    # the following covers the edge case where the old file was unversioned
    if is_new_version_of.version is None:
        is_new_version_of.version = previous_version
        is_new_version_of.save()
        if msg != "":
            msg += (
                f"& of new file to '{version}' (initial_version_id ="
                f" '{initial_version_id}')"
            )
    return new_file_id, initial_version_id, version  # type: ignore


def __init__(file: File, *args, **kwargs):
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
    data: Union[PathLike, DataLike] = kwargs.pop("data") if len(args) == 0 else args[0]
    key: Optional[str] = kwargs.pop("key") if "key" in kwargs else None
    run: Optional[Run] = kwargs.pop("run") if "run" in kwargs else None
    description: Optional[str] = (
        kwargs.pop("description") if "description" in kwargs else None
    )
    is_new_version_of: Optional[File] = (
        kwargs.pop("is_new_version_of") if "is_new_version_of" in kwargs else None
    )
    initial_version_id: Optional[str] = (
        kwargs.pop("initial_version_id") if "initial_version_id" in kwargs else None
    )
    version: Optional[str] = kwargs.pop("version") if "version" in kwargs else None
    format = kwargs.pop("format") if "format" in kwargs else None
    log_hint = kwargs.pop("log_hint") if "log_hint" in kwargs else True
    skip_check_exists = (
        kwargs.pop("skip_check_exists") if "skip_check_exists" in kwargs else False
    )

    if not len(kwargs) == 0:
        raise ValueError("Only data, key, run, description can be passed.")

    if is_new_version_of is None:
        provisional_id = init_id(version=version)
    else:
        if not isinstance(is_new_version_of, File):
            raise TypeError("is_new_version_of has to be of type ln.File")
        provisional_id, initial_version_id, version = get_ids_from_old_version_of_file(
            is_new_version_of, version
        )
        if description is None:
            description = is_new_version_of.description

    if version is not None:
        if initial_version_id is None:
            logger.info(
                "initializing versioning for this file! create future versions of it"
                " using ln.File(..., is_new_version_of=old_file)"
            )
    kwargs_or_file, privates = get_file_kwargs_from_data(
        data=data,
        key=key,
        run=run,
        format=format,
        provisional_id=provisional_id,
        skip_check_exists=skip_check_exists,
    )

    # an object with the same hash already exists
    if isinstance(kwargs_or_file, File):
        existing_file = kwargs_or_file
        # this is the way Django instantiates from the DB internally
        # https://github.com/django/django/blob/549d6ffeb6d626b023acc40c3bb2093b4b25b3d6/django/db/models/base.py#LL488C1-L491C51
        new_args = [
            getattr(existing_file, field.attname)
            for field in file._meta.concrete_fields
        ]
        super(File, file).__init__(*new_args)
        file._state.adding = False
        file._state.db = "default"
        return None
    else:
        kwargs = kwargs_or_file

    if isinstance(data, pd.DataFrame):
        if log_hint:
            logger.hint(
                "file is a dataframe, consider using File.from_df() to link column"
                " names as features"
            )
        kwargs["accessor"] = "DataFrame"
    elif data_is_anndata(data):
        if log_hint:
            logger.hint(
                "file is AnnDataLike, consider using File.from_anndata() to link"
                " var_names and obs.columns as features"
            )
        kwargs["accessor"] = "AnnData"
    elif data_is_mudata(data):
        kwargs["accessor"] = "MuData"

    kwargs["id"] = provisional_id
    kwargs["initial_version_id"] = initial_version_id
    kwargs["version"] = version
    kwargs["description"] = description
    # this check needs to come down here because key might be populated from an
    # existing file path during get_file_kwargs_from_data()
    if (
        kwargs["key"] is None
        and kwargs["description"] is None
        and kwargs["run"] is None
    ):
        raise ValueError("Pass one of key, run or description as a parameter")

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


@classmethod  # type: ignore
@doc_args(File.from_df.__doc__)
def from_df(
    cls,
    df: "pd.DataFrame",
    columns_ref: FieldAttr = Feature.name,
    key: Optional[str] = None,
    description: Optional[str] = None,
    run: Optional[Run] = None,
) -> "File":
    """{}"""
    file = File(data=df, key=key, run=run, description=description, log_hint=False)
    feature_set = FeatureSet.from_df(df)
    if feature_set is not None:
        file._feature_sets = {"columns": feature_set}
    else:
        file._feature_sets = {}
    return file


@classmethod  # type: ignore
@doc_args(File.from_anndata.__doc__)
def from_anndata(
    cls,
    adata: "AnnDataLike",
    var_ref: Optional[FieldAttr],
    obs_columns_ref: Optional[FieldAttr] = Feature.name,
    key: Optional[str] = None,
    description: Optional[str] = None,
    run: Optional[Run] = None,
) -> "File":
    """{}"""
    file = File(data=adata, key=key, run=run, description=description, log_hint=False)
    data_parse = adata
    if not isinstance(adata, AnnData):  # is a path
        filepath = adata if isinstance(adata, (Path, UPath)) else _str_to_path(adata)
        if isinstance(filepath, UPath):
            from lamindb.dev.storage._backed_access import backed_access

            data_parse = backed_access(filepath)
        else:
            data_parse = ad.read(filepath, backed="r")
        type = "float"
    else:
        type = convert_numpy_dtype_to_lamin_feature_type(adata.X.dtype)
    feature_sets = {}
    logger.info("parsing feature names of X stored in slot 'var'")
    logger.indent = "   "
    feature_set_var = FeatureSet.from_values(
        data_parse.var.index,
        var_ref,
        type=type,
    )

    if feature_set_var is not None:
        feature_sets["var"] = feature_set_var
        logger.save(f"linked: {feature_set_var}")
    logger.indent = ""
    if len(data_parse.obs.columns) > 0:
        logger.info("parsing feature names of slot 'obs'")
        logger.indent = "   "
        feature_set_obs = FeatureSet.from_df(data_parse.obs)
        if feature_set_obs is not None:
            feature_sets["obs"] = feature_set_obs
            logger.save(f"linked: {feature_set_obs}")
        logger.indent = ""
    file._feature_sets = feature_sets
    return file


@classmethod  # type: ignore
@doc_args(File.from_dir.__doc__)
def from_dir(
    cls,
    path: PathLike,
    key: Optional[str] = None,
    *,
    run: Optional[Run] = None,
) -> List["File"]:
    """{}"""
    folderpath = path if isinstance(path, (Path, UPath)) else _str_to_path(path)
    storage, use_existing_storage = process_pathlike(folderpath)
    # we are never erroring right now, can remove below, soon
    # if key is None and not use_existing_storage:
    #     formatted_roots = "\n".join(Storage.filter().list("root"))
    #     raise ValueError(
    #   "If `key` is None, don't support tracking folders outside one of the"
    #   f" storage roots:\n{formatted_roots}\nEither pass key, move folder"
    #   " or register new storage location:\nln.Storage(root='path/to/my/new_root',"
    #   " type='local').save()"
    #     )
    folder_key_path: Union[PurePath, Path]
    if key is None:
        if not use_existing_storage:
            logger.warning(
                "folder is outside existing storage location, will copy files from"
                f" {path} to {storage}/{folderpath.name}"
            )
            folder_key_path = Path(folderpath.name)
        else:
            # maintain the hierachy within an existing storage location
            folder_key_path = get_relative_path_to_directory(
                folderpath, storage.root_as_path()
            )
    else:
        folder_key_path = Path(key)

    # always sanitize by stripping a trailing slash
    folder_key = folder_key_path.as_posix().rstrip("/")

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
            # if creating from rglob, we don't need to check for existence
            file = File(filepath, run=run, key=file_key, skip_check_exists=True)
            files.append(file)
    settings.verbosity = verbosity
    logger.success(
        f"created {len(files)} files from directory using storage"
        f" {storage.root} and key = {folder_key}/"
    )
    return files


def replace(
    self,
    data: Union[PathLike, DataLike],
    run: Optional[Run] = None,
    format: Optional[str] = None,
) -> None:
    kwargs, privates = get_file_kwargs_from_data(
        provisional_id=self.id,
        data=data,
        key=self.key,
        run=run,
        format=format,
    )
    if self.key is not None:
        key_path = PurePosixPath(self.key)
        new_filename = f"{key_path.stem}{kwargs['suffix']}"
        # the following will only be true if the suffix changes!
        if key_path.name != new_filename:
            self._clear_storagekey = self.key
            self.key = str(key_path.with_name(new_filename))
            logger.warning(
                f"replacing the file will replace key '{key_path}' with '{self.key}'"
                f" and delete '{key_path}' upon `save()`"
            )
    else:
        self.key = kwargs["key"]
        old_storage = auto_storage_key_from_file(self)
        new_storage = (
            self.key if self.key is not None else f"{self.id}{kwargs['suffix']}"
        )
        if old_storage != new_storage:
            self._clear_storagekey = old_storage

    self.suffix = kwargs["suffix"]
    self.size = kwargs["size"]
    self.hash = kwargs["hash"]
    self.run = kwargs["run"]
    self._local_filepath = privates["local_filepath"]
    self._cloud_filepath = privates["cloud_filepath"]
    self._memory_rep = privates["memory_rep"]
    self._to_store = not privates[
        "check_path_in_storage"
    ]  # no need to upload if new file is already in storage


def backed(
    self, is_run_input: Optional[bool] = None
) -> Union["AnnDataAccessor", "BackedAccessor"]:
    suffixes = (".h5", ".hdf5", ".h5ad", ".zrad", ".zarr")
    if self.suffix not in suffixes:
        raise ValueError(
            "File should have a zarr or h5 object as the underlying data, please use"
            " one of the following suffixes for the object name:"
            f" {', '.join(suffixes)}."
        )

    from lamindb.dev.storage._backed_access import backed_access

    _track_run_input(self, is_run_input)

    filepath = filepath_from_file(self)
    # consider the case where an object is already locally cached
    localpath = setup_settings.instance.storage.cloud_to_local_no_update(filepath)
    if localpath.exists():
        return backed_access(localpath)
    else:
        return backed_access(filepath)


def _track_run_input(file: File, is_run_input: Optional[bool] = None):
    track_run_input = False
    if is_run_input is None:
        # we need a global run context for this to work
        if run_context.run is not None:
            # avoid cycles (a file is both input and output)
            if file.run != run_context.run:
                if settings.track_run_inputs:
                    transform_note = ""
                    if file.transform is not None:
                        transform_note = (
                            f", adding parent transform {file.transform.id}"
                        )
                    logger.info(
                        f"adding file {file.id} as input for run"
                        f" {run_context.run.id}{transform_note}"
                    )
                    track_run_input = True
                else:
                    logger.hint(
                        "track this file as a run input by passing `is_run_input=True`"
                    )
        else:
            if settings.track_run_inputs:
                logger.hint(
                    "you can auto-track this file as a run input by calling"
                    " `ln.track()`"
                )
    else:
        track_run_input = is_run_input
    if track_run_input:
        if run_context.run is None:
            raise ValueError(
                "No global run context set. Call ln.context.track() or link input to a"
                " run object via `run.input_files.append(file)`"
            )
        # avoid adding the same run twice
        # avoid cycles (a file is both input and output)
        if not file.input_of.contains(run_context.run) and file.run != run_context.run:
            run_context.run.save()
            file.input_of.add(run_context.run)
            run_context.run.transform.parents.add(file.transform)


def load(self, is_run_input: Optional[bool] = None, stream: bool = False) -> DataLike:
    _track_run_input(self, is_run_input)
    if hasattr(self, "_memory_rep") and self._memory_rep is not None:
        return self._memory_rep
    return load_to_memory(filepath_from_file(self), stream=stream)


def stage(self, is_run_input: Optional[bool] = None) -> Path:
    if self.suffix in (".zrad", ".zarr"):
        raise RuntimeError("zarr object can't be staged, please use load() or stream()")
    _track_run_input(self, is_run_input)

    filepath = filepath_from_file(self)
    cb = ProgressCallback("downloading")
    return setup_settings.instance.storage.cloud_to_local(filepath, callback=cb)


def delete(self, storage: Optional[bool] = None) -> None:
    if storage is None:
        response = input(f"Are you sure you want to delete {self} from storage? (y/n)")
        delete_in_storage = response == "y"
    else:
        delete_in_storage = storage

    # need to grab file path before deletion
    filepath = self.path
    # only delete in storage if DB delete is successful
    # DB delete might error because of a foreign key constraint violated etc.
    self._delete_skip_storage()
    # we don't yet have any way to bring back the deleted metadata record
    # in case the storage deletion fails - this is important for ACID down the road
    if delete_in_storage:
        delete_storage(filepath)
        logger.success(f"deleted stored object {colors.yellow(f'{filepath}')}")


def _delete_skip_storage(file, *args, **kwargs) -> None:
    super(File, file).delete(*args, **kwargs)


def save(self, *args, **kwargs) -> None:
    self._save_skip_storage(*args, **kwargs)
    from lamindb._save import check_and_attempt_clearing, check_and_attempt_upload

    exception = check_and_attempt_upload(self)
    if exception is not None:
        self._delete_skip_storage()
        raise RuntimeError(exception)
    exception = check_and_attempt_clearing(self)
    if exception is not None:
        raise RuntimeError(exception)


def _save_skip_storage(file, *args, **kwargs) -> None:
    if file.transform is not None:
        file.transform.save()
    if file.run is not None:
        file.run.save()
    if hasattr(file, "_feature_sets"):
        for feature_set in file._feature_sets.values():
            feature_set.save()
        s = "s" if len(file._feature_sets) > 1 else ""
        logger.save(
            f"saved {len(file._feature_sets)} feature set{s} for slot{s}:"
            f" {list(file._feature_sets.keys())}"
        )
    super(File, file).save(*args, **kwargs)
    if hasattr(file, "_feature_sets"):
        links = []
        for slot, feature_set in file._feature_sets.items():
            links.append(
                File.feature_sets.through(
                    file_id=file.id, feature_set_id=feature_set.id, slot=slot
                )
            )

        from lamindb._save import bulk_create

        bulk_create(links)


@property  # type: ignore
@doc_args(File.path.__doc__)
def path(self) -> Union[Path, UPath]:
    """{}"""
    return filepath_from_file(self)


# adapted from: https://stackoverflow.com/questions/9727673
@classmethod  # type: ignore
@doc_args(File.tree.__doc__)
def tree(
    cls: File,
    path: Optional[PathLike] = None,
    *,
    level: int = -1,
    limit_to_directories: bool = False,
    length_limit: int = 1000,
) -> None:
    """{}"""
    space = "    "
    branch = "│   "
    tee = "├── "
    last = "└── "

    if path is None:
        dir_path = settings.storage
    else:
        dir_path = path if isinstance(path, (Path, UPath)) else _str_to_path(path)
    n_files = 0
    n_directories = 0

    # by default only including registered files
    # need a flag and a proper implementation
    registered_paths: Set[Any] = set()
    registered_dirs: Set[Any] = set()
    if path is None:
        registered_paths = {
            file.path for file in cls.filter(storage_id=setup_settings.storage.id).all()
        }
        registered_dirs = {d for p in registered_paths for d in p.parents}

    def inner(dir_path: Union[Path, UPath], prefix: str = "", level=-1):
        nonlocal n_files, n_directories
        if not level:
            return  # 0, stop iterating
        stripped_dir_path = dir_path.as_posix().rstrip("/")
        # do not iterate through zarr directories
        if stripped_dir_path.endswith((".zarr", ".zrad")):
            return
        # this is needed so that the passed folder is not listed
        contents = [
            i
            for i in dir_path.iterdir()
            if i.as_posix().rstrip("/") != stripped_dir_path
        ]
        if limit_to_directories:
            contents = [d for d in contents if d.is_dir()]
        pointers = [tee] * (len(contents) - 1) + [last]
        for pointer, path in zip(pointers, contents):
            if path.is_dir():
                if registered_dirs and path not in registered_dirs:
                    continue
                yield prefix + pointer + path.name
                n_directories += 1
                extension = branch if pointer == tee else space
                yield from inner(path, prefix=prefix + extension, level=level - 1)
            elif not limit_to_directories:
                if registered_paths and path not in registered_paths:
                    continue
                yield prefix + pointer + path.name
                n_files += 1

    folder_tree = ""
    iterator = inner(dir_path, level=level)
    for line in islice(iterator, length_limit):
        folder_tree += f"\n{line}"
    if next(iterator, None):
        folder_tree += f"... length_limit, {length_limit}, reached, counted:"
    directory_info = "directory" if n_directories == 1 else "directories"
    print(
        f"{dir_path.name} ({n_directories} sub-{directory_info} & {n_files} files):"
        f" {folder_tree}"
    )


def inherit_relations(self, file: File, fields: Optional[List[str]] = None):
    """Inherit many-to-many relationships from another file.

    Examples:
        >>> file1 = ln.File(pd.DataFrame(index=[0,1]))
        >>> file1.save()
        >>> file2 = ln.File(pd.DataFrame(index=[2,3]))
        >>> file2.save()
        >>> ln.save(ln.Label.from_values(["Label1", "Label2", "Label3"], field="name"))
        >>> labels = ln.Label.filter(name__icontains = "label").all()
        >>> file1.labels.set(labels)
        >>> file2.inherit_relations(file1, ["labels"])
        💬 Inheriting 1 field: ['labels']
        >>> file2.labels.list("name")
        ['Label1', 'Label2', 'Label3']
    """
    if fields is None:
        # fields in the model definition
        related_names = [i.name for i in file._meta.many_to_many]
        # fields back linked
        related_names += [i.related_name for i in file._meta.related_objects]
    else:
        related_names = []
        for field in fields:
            if hasattr(file, field):
                related_names.append(field)
            else:
                raise KeyError(f"No many-to-many relationship is found with '{field}'")

    if None in related_names:
        related_names.remove(None)

    inherit_names = [
        related_name
        for related_name in related_names
        if related_name is not None
        if file.__getattribute__(related_name).exists()
    ]

    s = "s" if len(inherit_names) > 1 else ""
    logger.info(f"inheriting {len(inherit_names)} field{s}: {inherit_names}")
    for related_name in inherit_names:
        self.__getattribute__(related_name).set(
            file.__getattribute__(related_name).all()
        )


METHOD_NAMES = [
    "__init__",
    "from_anndata",
    "from_df",
    "backed",
    "stage",
    "load",
    "delete",
    "save",
    "replace",
    "from_dir",
    "tree",
]

if _TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(File, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, File, globals())

# privates currently dealt with separately
File._delete_skip_storage = _delete_skip_storage
File._save_skip_storage = _save_skip_storage
# TODO: move these to METHOD_NAMES
setattr(File, "view_lineage", view_lineage)
setattr(File, "inherit_relations", inherit_relations)
setattr(File, "path", path)
