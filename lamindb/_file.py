from pathlib import Path, PurePath, PurePosixPath
from typing import Any, List, Optional, Tuple, Union

import anndata as ad
import fsspec
import lamindb_setup
import pandas as pd
from anndata import AnnData
from lamin_utils import colors, logger
from lamindb_setup import settings as setup_settings
from lamindb_setup._init_instance import register_storage
from lamindb_setup.dev import StorageSettings
from lamindb_setup.dev._docs import doc_args
from lamindb_setup.dev._hub_utils import get_storage_region
from lamindb_setup.dev.upath import create_path, extract_suffix_from_path
from lnschema_core import Feature, FeatureSet, File, Run, Storage
from lnschema_core.models import IsTree
from lnschema_core.types import (
    AnnDataLike,
    DataLike,
    FieldAttr,
    PathLike,
    VisibilityChoice,
)

from lamindb._utils import attach_func_to_class_method
from lamindb.dev._data import _track_run_input
from lamindb.dev._settings import settings
from lamindb.dev.hashing import b16_to_b64, hash_file
from lamindb.dev.storage import (
    LocalPathClasses,
    UPath,
    delete_storage,
    infer_suffix,
    load_to_memory,
    size_adata,
    write_to_file,
)
from lamindb.dev.storage._backed_access import AnnDataAccessor, BackedAccessor
from lamindb.dev.storage.file import (
    auto_storage_key_from_file,
    auto_storage_key_from_id_suffix,
    filepath_from_file,
)
from lamindb.dev.versioning import get_ids_from_old_version, init_uid

from . import _TESTING
from ._feature import convert_numpy_dtype_to_lamin_feature_type
from .dev._data import (
    add_transform_to_kwargs,
    get_run,
    save_feature_set_links,
    save_feature_sets,
)
from .dev.storage.file import AUTO_KEY_PREFIX


def process_pathlike(
    filepath: UPath, skip_existence_check: bool = False
) -> Tuple[Storage, bool]:
    if not skip_existence_check:
        try:  # check if file exists
            if not filepath.exists():
                raise FileNotFoundError(filepath)
        except PermissionError:
            pass
    if isinstance(filepath, LocalPathClasses):
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
            if not isinstance(filepath, LocalPathClasses):
                # for a cloud path, new_root is always the bucket name
                new_root = list(filepath.parents)[-1]
                new_root_str = new_root.as_posix().rstrip("/")
                region = get_storage_region(new_root_str)
                storage_settings = StorageSettings(new_root_str, region)
                storage_record = register_storage(storage_settings)
                use_existing_storage_key = True
                return storage_record, use_existing_storage_key
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
    provisional_uid: str,
    data: Union[PathLike, DataLike],
    format: Optional[str],
    key: Optional[str],
    skip_existence_check: bool = False,
) -> Tuple[Any, Union[Path, UPath], str, Storage, bool]:
    """Serialize a data object that's provided as file or in memory."""
    # if not overwritten, data gets stored in default storage
    if isinstance(data, (str, Path, UPath)):  # PathLike, spelled out
        filepath = create_path(data)
        storage, use_existing_storage_key = process_pathlike(
            filepath, skip_existence_check=skip_existence_check
        )
        suffix = extract_suffix_from_path(filepath)
        memory_rep = None
    elif isinstance(data, (pd.DataFrame, AnnData)):  # DataLike, spelled out
        storage = lamindb_setup.settings.storage.record
        memory_rep = data
        if key is not None:
            key_suffix = extract_suffix_from_path(PurePosixPath(key), arg_name="key")
            # use suffix as the (adata) format if the format is not provided
            if isinstance(data, AnnData) and format is None and len(key_suffix) > 0:
                format = key_suffix[1:]
        else:
            key_suffix = None
        suffix = infer_suffix(data, format)
        if key_suffix is not None and key_suffix != suffix:
            raise ValueError(
                f"The suffix '{key_suffix}' of the provided key is incorrect, it should"
                f" be '{suffix}'."
            )
        cache_name = f"{provisional_uid}{suffix}"
        filepath = lamindb_setup.settings.storage.cache_dir / cache_name
        # Alex: I don't understand the line below
        if filepath.suffixes == []:
            filepath = filepath.with_suffix(suffix)
        if suffix not in {".zarr", ".zrad"}:
            write_to_file(data, filepath)
        use_existing_storage_key = False
    else:
        raise NotImplementedError(
            f"Do not know how to create a file object from {data}, pass a filepath"
            " instead!"
        )
    return memory_rep, filepath, suffix, storage, use_existing_storage_key


def get_hash(
    filepath: UPath,
    suffix,
    filepath_stat=None,
    check_hash: bool = True,
) -> Union[Tuple[Optional[str], Optional[str]], File]:
    if suffix in {".zarr", ".zrad"}:
        return None
    if not isinstance(filepath, LocalPathClasses):
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
    # also checks hidden and trashed files
    result = File.filter(hash=hash, visibility=None).list()
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
            if result[0].visibility < 1:
                if result[0].visibility == -1:
                    visibility_text = "in the trash"
                elif result[0].visibility == 0:
                    visibility_text = "hidden"
                logger.warning(
                    f"the existing file is {visibility_text}, restore it before use:"
                    " `file.restore()`"
                )
            return result[0]
    else:
        return hash, hash_type


def get_path_size_hash(
    filepath: UPath,
    memory_rep: Optional[Union[pd.DataFrame, AnnData]],
    suffix: str,
    check_hash: bool = True,
):
    cloudpath = None
    localpath = None
    hash_and_type: Tuple[Optional[str], Optional[str]]

    if suffix in {".zarr", ".zrad"}:
        if memory_rep is not None:
            size = size_adata(memory_rep)
        else:
            if not isinstance(filepath, LocalPathClasses):
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
            if not isinstance(filepath, LocalPathClasses):
                size = filepath_stat["size"]
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
        if check_path_is_child_of_root(filepath, root=create_path(storage.root)):
            return storage
    return False


def check_path_is_child_of_root(
    filepath: Union[Path, UPath], root: Optional[Union[Path, UPath]] = None
) -> bool:
    if root is None:
        root = lamindb_setup.settings.storage.root

    filepath = UPath(str(filepath)) if not isinstance(filepath, UPath) else filepath
    root = UPath(str(root)) if not isinstance(root, UPath) else root

    # the following comparisons can fail if types aren't comparable
    if not isinstance(filepath, LocalPathClasses) and not isinstance(
        root, LocalPathClasses
    ):
        # the following tests equivalency of two UPath objects
        # via string representations; otherwise
        # S3Path('s3://lndb-storage/') and S3Path('s3://lamindb-ci/')
        # test as equivalent
        return list(filepath.parents)[-1].as_posix() == root.as_posix()
    elif isinstance(filepath, LocalPathClasses) and isinstance(root, LocalPathClasses):
        return root.resolve() in filepath.resolve().parents
    else:
        return False


def get_relative_path_to_directory(
    path: Union[PurePath, Path, UPath], directory: Union[PurePath, Path, UPath]
) -> Union[PurePath, Path]:
    if isinstance(directory, UPath) and not isinstance(directory, LocalPathClasses):
        # UPath.relative_to() is not behaving as it should (2023-04-07)
        # need to lstrip otherwise inconsistent behavior across trailing slashes
        # see test_file.py: test_get_relative_path_to_directory
        relpath = PurePath(
            path.as_posix().replace(directory.as_posix(), "").lstrip("/")
        )
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
    provisional_uid: str,
    skip_check_exists: bool = False,
):
    run = get_run(run)
    memory_rep, filepath, suffix, storage, use_existing_storage_key = process_data(
        provisional_uid, data, format, key, skip_check_exists
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
    if use_existing_storage_key:
        inferred_key = get_relative_path_to_directory(
            path=filepath, directory=storage.root_as_path()
        ).as_posix()
        if key is None:
            key = inferred_key
        else:
            if not key == inferred_key:
                raise ValueError(
                    f"The path '{data}' is already in registered storage"
                    f" '{storage.root}' with key '{inferred_key}'\nYou passed"
                    f" conflicting key '{key}': please move the file before"
                    " registering it."
                )
        check_path_in_storage = True
    else:
        storage = lamindb_setup.settings.storage.record

    if key is not None and key.startswith(AUTO_KEY_PREFIX):
        raise ValueError(f"Key cannot start with {AUTO_KEY_PREFIX}")

    log_storage_hint(
        check_path_in_storage=check_path_in_storage,
        storage=storage,
        key=key,
        uid=provisional_uid,
        suffix=suffix,
    )

    # do we use a virtual or an actual storage key?
    key_is_virtual = settings.file_use_virtual_keys

    # if the file is already in storage, independent of the default
    # we use an actual storage key
    if check_path_in_storage:
        key_is_virtual = False

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
        key_is_virtual=key_is_virtual,
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
    uid: str,
    suffix: str,
) -> None:
    hint = ""
    if check_path_in_storage:
        display_root = storage.root  # type: ignore
        # check whether path is local
        if fsspec.utils.get_protocol(storage.root) == "file":  # type: ignore
            # if it's a local path, check whether it's in the current working directory
            root_path = Path(storage.root)  # type: ignore
            if check_path_is_child_of_root(root_path, Path.cwd()):
                # only display the relative path, not the fully resolved path
                display_root = root_path.relative_to(Path.cwd())
        hint += f"file in storage '{display_root}'"  # type: ignore
    else:
        hint += "file will be copied to default storage upon `save()`"
    if key is None:
        storage_key = auto_storage_key_from_id_suffix(uid, suffix)
        hint += f" with key `None` ('{storage_key}')"
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
    initial_version_id: Optional[int] = (
        kwargs.pop("initial_version_id") if "initial_version_id" in kwargs else None
    )
    version: Optional[str] = kwargs.pop("version") if "version" in kwargs else None
    visibility: Optional[int] = (
        kwargs.pop("visibility")
        if "visibility" in kwargs
        else VisibilityChoice.default.value
    )
    format = kwargs.pop("format") if "format" in kwargs else None
    log_hint = kwargs.pop("log_hint") if "log_hint" in kwargs else True
    skip_check_exists = (
        kwargs.pop("skip_check_exists") if "skip_check_exists" in kwargs else False
    )

    if not len(kwargs) == 0:
        raise ValueError(
            "Only data, key, run, description, version, is_new_version_of, visibility"
            f" can be passed, you passed: {kwargs}"
        )

    if is_new_version_of is None:
        provisional_uid = init_uid(version=version, n_full_id=20)
    else:
        if not isinstance(is_new_version_of, File):
            raise TypeError("is_new_version_of has to be of type ln.File")
        provisional_uid, initial_version_id, version = get_ids_from_old_version(
            is_new_version_of, version, n_full_id=20
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
        provisional_uid=provisional_uid,
        skip_check_exists=skip_check_exists,
    )

    # an object with the same hash already exists
    if isinstance(kwargs_or_file, File):
        from ._registry import init_self_from_db

        # kwargs_or_file is an existing file
        init_self_from_db(file, kwargs_or_file)
        return None
    else:
        kwargs = kwargs_or_file

    if isinstance(data, pd.DataFrame):
        if log_hint:
            logger.hint(
                "data is a dataframe, consider using .from_df() to link column"
                " names as features"
            )
        kwargs["accessor"] = "DataFrame"
    elif data_is_anndata(data):
        if log_hint:
            logger.hint(
                "data is AnnDataLike, consider using .from_anndata() to link"
                " var_names and obs.columns as features"
            )
        kwargs["accessor"] = "AnnData"
    elif data_is_mudata(data):
        kwargs["accessor"] = "MuData"

    kwargs["uid"] = provisional_uid
    kwargs["initial_version_id"] = initial_version_id
    kwargs["version"] = version
    kwargs["description"] = description
    kwargs["visibility"] = visibility
    # this check needs to come down here because key might be populated from an
    # existing file path during get_file_kwargs_from_data()
    if (
        kwargs["key"] is None
        and kwargs["description"] is None
        and kwargs["run"] is None
    ):
        raise ValueError("Pass one of key, run or description as a parameter")

    add_transform_to_kwargs(kwargs, kwargs["run"])

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
    field: FieldAttr = Feature.name,
    key: Optional[str] = None,
    description: Optional[str] = None,
    run: Optional[Run] = None,
    version: Optional[str] = None,
    is_new_version_of: Optional["File"] = None,
    **kwargs,
) -> "File":
    """{}"""
    file = File(
        data=df,
        key=key,
        run=run,
        description=description,
        version=version,
        is_new_version_of=is_new_version_of,
        log_hint=False,
    )
    feature_set = FeatureSet.from_df(df, field=field, **kwargs)
    if feature_set is not None:
        file._feature_sets = {"columns": feature_set}
    else:
        file._feature_sets = {}
    return file


def parse_feature_sets_from_anndata(
    adata: AnnDataLike,
    field: Optional[FieldAttr],
    **kwargs,
):
    data_parse = adata
    if not isinstance(adata, AnnData):  # is a path
        filepath = create_path(adata)  # returns Path for local
        if not isinstance(filepath, LocalPathClasses):
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
        field,
        type=type,
        **kwargs,
    )
    if feature_set_var is not None:
        feature_sets["var"] = feature_set_var
        logger.save(f"linked: {feature_set_var}")
    logger.indent = ""
    if len(data_parse.obs.columns) > 0:
        logger.info("parsing feature names of slot 'obs'")
        logger.indent = "   "
        feature_set_obs = FeatureSet.from_df(
            data_parse.obs,
            **kwargs,
        )
        if feature_set_obs is not None:
            feature_sets["obs"] = feature_set_obs
            logger.save(f"linked: {feature_set_obs}")
        logger.indent = ""
    return feature_sets


@classmethod  # type: ignore
@doc_args(File.from_anndata.__doc__)
def from_anndata(
    cls,
    adata: "AnnDataLike",
    field: Optional[FieldAttr],
    key: Optional[str] = None,
    description: Optional[str] = None,
    run: Optional[Run] = None,
    version: Optional[str] = None,
    is_new_version_of: Optional["File"] = None,
    **kwargs,
) -> "File":
    """{}"""
    file = File(
        data=adata,
        key=key,
        run=run,
        description=description,
        version=version,
        is_new_version_of=is_new_version_of,
        log_hint=False,
    )
    file._feature_sets = parse_feature_sets_from_anndata(adata, field, **kwargs)
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
    folderpath: UPath = create_path(path)  # returns Path for local
    storage, use_existing_storage = process_pathlike(folderpath)
    folder_key_path: Union[PurePath, Path]
    if key is None:
        if not use_existing_storage:
            logger.warning(
                "folder is outside existing storage location, will copy files from"
                f" {path} to {storage.root}/{folderpath.name}"
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

    # TODO: (non-local) UPath doesn't list the first level files and dirs with "*"
    pattern = "" if not isinstance(folderpath, LocalPathClasses) else "*"

    # silence fine-grained logging
    verbosity = settings.verbosity
    verbosity_int = settings._verbosity_int
    if verbosity_int >= 1:
        settings.verbosity = "warning"
    files_dict = {}
    for filepath in folderpath.rglob(pattern):
        if filepath.is_file():
            relative_path = get_relative_path_to_directory(filepath, folderpath)
            file_key = folder_key + "/" + relative_path.as_posix()
            # if creating from rglob, we don't need to check for existence
            file = File(filepath, run=run, key=file_key, skip_check_exists=True)
            files_dict[file.uid] = file
    settings.verbosity = verbosity

    # run sanity check on hashes
    hashes = [file.hash for file in files_dict.values() if file.hash is not None]
    uids = files_dict.keys()
    if len(set(hashes)) == len(hashes):
        files = list(files_dict.values())
    else:
        # consider exact duplicates (same id, same hash)
        # below can't happen anymore because files is a dict now
        # if len(set(uids)) == len(set(hashes)):
        #     logger.warning("dropping duplicate records in list of file records")
        #     files = list(set(uids))
        # consider false duplicates (different id, same hash)
        if not len(set(uids)) == len(set(hashes)):
            seen_hashes = set()
            non_unique_files = {
                hash: file
                for hash, file in files_dict.items()
                if file.hash in seen_hashes or seen_hashes.add(file.hash)  # type: ignore  # noqa
            }
            display_non_unique = "\n    ".join(f"{file}" for file in non_unique_files)
            logger.warning(
                "there are multiple file uids with the same hashes, dropping"
                f" {len(non_unique_files)} duplicates out of {len(files_dict)} files:\n"
                f"    {display_non_unique}"
            )
            files = [
                file
                for file in files_dict.values()
                if file not in non_unique_files.values()
            ]
    logger.success(
        f"created {len(files)} files from directory using storage"
        f" {storage.root} and key = {folder_key}/"
    )
    return files


# docstring handled through attach_func_to_class_method
def replace(
    self,
    data: Union[PathLike, DataLike],
    run: Optional[Run] = None,
    format: Optional[str] = None,
) -> None:
    kwargs, privates = get_file_kwargs_from_data(
        provisional_uid=self.uid,
        data=data,
        key=self.key,
        run=run,
        format=format,
    )

    # this file already exists
    if privates is None:
        return kwargs

    check_path_in_storage = privates["check_path_in_storage"]
    if check_path_in_storage:
        raise ValueError("Can only replace with a local file not in any Storage.")

    if self.key is not None and not self.key_is_virtual:
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
        old_storage = auto_storage_key_from_file(self)
        new_storage = auto_storage_key_from_id_suffix(self.uid, kwargs["suffix"])
        if old_storage != new_storage:
            self._clear_storagekey = old_storage
            if self.key is not None:
                new_key_path = PurePosixPath(self.key).with_suffix(kwargs["suffix"])
                self.key = str(new_key_path)

    self.suffix = kwargs["suffix"]
    self.size = kwargs["size"]
    self.hash = kwargs["hash"]
    self.hash_type = kwargs["hash_type"]
    self.run_id = kwargs["run_id"]
    self.run = kwargs["run"]

    self._local_filepath = privates["local_filepath"]
    self._cloud_filepath = privates["cloud_filepath"]
    self._memory_rep = privates["memory_rep"]
    # no need to upload if new file is already in storage
    self._to_store = not check_path_in_storage


# docstring handled through attach_func_to_class_method
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


# docstring handled through attach_func_to_class_method
def load(
    self, is_run_input: Optional[bool] = None, stream: bool = False, **kwargs
) -> DataLike:
    _track_run_input(self, is_run_input)
    if hasattr(self, "_memory_rep") and self._memory_rep is not None:
        return self._memory_rep
    return load_to_memory(filepath_from_file(self), stream=stream, **kwargs)


# docstring handled through attach_func_to_class_method
def stage(self, is_run_input: Optional[bool] = None) -> Path:
    if self.suffix in {".zrad", ".zarr"}:
        raise RuntimeError("zarr object can't be staged, please use load() or stream()")
    _track_run_input(self, is_run_input)

    filepath = filepath_from_file(self)
    return setup_settings.instance.storage.cloud_to_local(filepath, print_progress=True)


# docstring handled through attach_func_to_class_method
def delete(
    self, permanent: Optional[bool] = None, storage: Optional[bool] = None
) -> None:
    # by default, we only move files into the trash
    if self.visibility > VisibilityChoice.trash.value and permanent is not True:
        if storage is not None:
            logger.warning("moving file to trash, storage arg is ignored")
        # move to trash
        self.visibility = VisibilityChoice.trash.value
        self.save()
        logger.warning("moved file to trash")
        return

    # if the file is already in the trash
    # permanent delete skips the trash
    if permanent is None:
        response = input(
            "File record is already in trash! Are you sure you want to permanently"
            " delete it? (y/n) You can't undo this action."
        )
        delete_record = response == "y"
    else:
        delete_record = permanent

    if delete_record:
        # need to grab file path before deletion
        filepath = self.path
        # only delete in storage if DB delete is successful
        # DB delete might error because of a foreign key constraint violated etc.
        self._delete_skip_storage()
        if self.key is None or self.key_is_virtual:
            delete_in_storage = True
            if storage is not None:
                logger.warning("storage arg is ignored if storage key is non-semantic")
        else:
            # for files with non-virtual semantic storage keys (key is not None)
            # ask for extra-confirmation
            if storage is None:
                response = input(
                    f"Are you sure to want to delete {filepath}? (y/n)  You can't undo"
                    " this action."
                )
                delete_in_storage = response == "y"
            else:
                delete_in_storage = storage
        # we don't yet have logic to bring back the deleted metadata record
        # in case storage deletion fails - this is important for ACID down the road
        if delete_in_storage:
            delete_storage(filepath)
            logger.success(f"deleted {colors.yellow(f'{filepath}')}")


def _delete_skip_storage(file, *args, **kwargs) -> None:
    super(File, file).delete(*args, **kwargs)


# docstring handled through attach_func_to_class_method
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
    save_feature_sets(file)
    super(File, file).save(*args, **kwargs)
    save_feature_set_links(file)


@property  # type: ignore
@doc_args(File.path.__doc__)
def path(self) -> Union[Path, UPath]:
    """{}"""
    return filepath_from_file(self)


@classmethod  # type: ignore
@doc_args(IsTree.view_tree.__doc__)
def view_tree(
    cls,
    level: int = -1,
    limit_to_directories: bool = False,
    length_limit: int = 1000,
    max_files_per_dir_per_type: int = 7,
) -> None:
    """{}"""
    from lamindb.dev._view_tree import view_tree as _view_tree

    _view_tree(
        cls=cls,
        level=level,
        limit_to_directories=limit_to_directories,
        length_limit=length_limit,
        max_files_per_dir_per_type=max_files_per_dir_per_type,
    )


# docstring handled through attach_func_to_class_method
def restore(self) -> None:
    self.visibility = VisibilityChoice.default.value
    self.save()


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
    "restore",
    "view_tree",
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
setattr(File, "path", path)
# this seems a Django-generated function
delattr(File, "get_visibility_display")
