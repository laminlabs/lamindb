from __future__ import annotations

import shutil
from pathlib import Path, PurePath, PurePosixPath
from typing import TYPE_CHECKING, Any, Mapping

import fsspec
import lamindb_setup as ln_setup
import pandas as pd
from anndata import AnnData
from lamin_utils import colors, logger
from lamindb_setup import settings as setup_settings
from lamindb_setup._init_instance import register_storage_in_instance
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core._settings_storage import init_storage
from lamindb_setup.core.hashing import hash_dir, hash_file
from lamindb_setup.core.upath import (
    create_path,
    extract_suffix_from_path,
    get_stat_dir_cloud,
    get_stat_file_cloud,
)
from lnschema_core.models import Artifact, FeatureManager, ParamManager, Run, Storage
from lnschema_core.types import (
    VisibilityChoice,
)

from lamindb._utils import attach_func_to_class_method
from lamindb.core._data import HasFeatures, _track_run_input
from lamindb.core._settings import settings
from lamindb.core.storage import (
    LocalPathClasses,
    UPath,
    delete_storage,
    infer_suffix,
    load_to_memory,
    write_to_disk,
)
from lamindb.core.storage.paths import (
    auto_storage_key_from_artifact,
    auto_storage_key_from_artifact_uid,
    filepath_from_artifact,
)
from lamindb.core.versioning import get_uid_from_old_version, init_uid

from .core._data import (
    add_transform_to_kwargs,
    get_run,
    save_feature_set_links,
    save_feature_sets,
)
from .core.storage.objects import _mudata_is_installed
from .core.storage.paths import AUTO_KEY_PREFIX

try:
    from .core.storage._zarr import zarr_is_adata
except ImportError:

    def zarr_is_adata(storepath):  # type: ignore
        raise ImportError("Please install zarr: pip install zarr")


if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr
    from mudata import MuData
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment

    from lamindb.core.storage._backed_access import AnnDataAccessor, BackedAccessor


def process_pathlike(
    filepath: UPath,
    default_storage: Storage,
    using_key: str | None,
    skip_existence_check: bool = False,
) -> tuple[Storage, bool]:
    if not skip_existence_check:
        try:  # check if file exists
            if not filepath.exists():
                raise FileNotFoundError(filepath)
        except PermissionError:
            pass
    if isinstance(filepath, LocalPathClasses):
        filepath = filepath.resolve()
    if check_path_is_child_of_root(filepath, default_storage.root):
        use_existing_storage_key = True
        return default_storage, use_existing_storage_key
    else:
        # check whether the path is part of one of the existing
        # already-registered storage locations
        result = False
        # within the hub, we don't want to perform check_path_in_existing_storage
        if using_key is None:
            result = check_path_in_existing_storage(filepath, using_key)
        if isinstance(result, Storage):
            use_existing_storage_key = True
            return result, use_existing_storage_key
        else:
            # if the path is in the cloud, we have a good candidate
            # for the storage root: the bucket
            if not isinstance(filepath, LocalPathClasses):
                # for a cloud path, new_root is always the bucket name
                new_root = list(filepath.parents)[-1]
                storage_settings = init_storage(new_root)
                storage_record = register_storage_in_instance(storage_settings)
                use_existing_storage_key = True
                return storage_record, use_existing_storage_key
            # if the filepath is local
            else:
                use_existing_storage_key = False
                # if the default storage is local we'll throw an error if the user
                # doesn't provide a key
                if default_storage.type == "local":
                    return default_storage, use_existing_storage_key
                # if the default storage is in the cloud (the file is going to
                # be uploaded upon saving it), we treat the filepath as a cache
                else:
                    return default_storage, use_existing_storage_key


def process_data(
    provisional_uid: str,
    data: UPathStr | pd.DataFrame | AnnData,
    format: str | None,
    key: str | None,
    default_storage: Storage,
    using_key: str | None,
    skip_existence_check: bool = False,
) -> tuple[Any, Path | UPath, str, Storage, bool]:
    """Serialize a data object that's provided as file or in memory."""
    # if not overwritten, data gets stored in default storage
    if _mudata_is_installed():
        from mudata import MuData

        data_types = (pd.DataFrame, AnnData, MuData)
    else:
        data_types = (pd.DataFrame, AnnData)  # type:ignore

    if isinstance(data, (str, Path, UPath)):  # UPathStr, spelled out
        access_token = (
            default_storage._access_token
            if hasattr(default_storage, "_access_token")
            else None
        )
        path = create_path(data, access_token=access_token)
        storage, use_existing_storage_key = process_pathlike(
            path,
            default_storage=default_storage,
            using_key=using_key,
            skip_existence_check=skip_existence_check,
        )
        suffix = extract_suffix_from_path(path)
        memory_rep = None
    elif isinstance(data, data_types):
        storage = default_storage
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
        path = settings._storage_settings.cache_dir / cache_name
        # Alex: I don't understand the line below
        if path.suffixes == []:
            path = path.with_suffix(suffix)
        write_to_disk(data, path)
        use_existing_storage_key = False
    else:
        raise NotImplementedError(
            f"Do not know how to create a artifact object from {data}, pass a path"
            " instead!"
        )
    return memory_rep, path, suffix, storage, use_existing_storage_key


def get_stat_or_artifact(
    path: UPath,
    check_hash: bool = True,
    using_key: str | None = None,
) -> tuple[int, str | None, str | None, int | None] | Artifact:
    n_objects = None
    if settings.creation.artifact_skip_size_hash:
        return None, None, None, n_objects
    stat = path.stat()  # one network request
    if not isinstance(path, LocalPathClasses):
        size, hash, hash_type = None, None, None
        if stat is not None:
            # convert UPathStatResult to fsspec info dict
            stat = stat.as_info()
            if "ETag" in stat:  # is file
                size, hash, hash_type = get_stat_file_cloud(stat)
            elif stat["type"] == "directory":
                size, hash, hash_type, n_objects = get_stat_dir_cloud(path)
        if hash is None:
            logger.warning(f"did not add hash for {path}")
            return size, hash, hash_type, n_objects
    else:
        if path.is_dir():
            size, hash, hash_type, n_objects = hash_dir(path)
        else:
            hash, hash_type = hash_file(path)
            size = stat.st_size
    if not check_hash:
        return size, hash, hash_type, n_objects
    # also checks hidden and trashed files
    # in Alex's mind the following two lines should be equivalent
    # but they aren't according to pytest tests/test_artifact.py::test_from_dir_single_artifact
    if using_key is None:
        result = Artifact.filter(hash=hash, visibility=None).all()
    else:
        result = (
            Artifact.objects.using(using_key).filter(hash=hash, visibility=None).all()
        )
    if len(result) > 0:
        if settings.creation.artifact_if_hash_exists == "error":
            msg = f"artifact with same hash exists: {result[0]}"
            hint = (
                "💡 you can make this error a warning:\n"
                "    ln.settings.creation.artifact_if_hash_exists"
            )
            raise FileExistsError(f"{msg}\n{hint}")
        elif settings.creation.artifact_if_hash_exists == "warn_create_new":
            logger.warning(
                "creating new Artifact object despite existing artifact with same hash:"
                f" {result[0]}"
            )
            return size, hash, hash_type, n_objects
        else:
            if result[0].visibility == -1:
                raise FileExistsError(
                    f"You're trying to re-create this artifact in trash: {result[0]}"
                    "Either permanently delete it with `artifact.delete(permanent=True)` or restore it with `artifact.restore()`"
                )
            logger.important(f"returning existing artifact with same hash: {result[0]}")
            return result[0]
    else:
        return size, hash, hash_type, n_objects


def check_path_in_existing_storage(
    path: Path | UPath, using_key: str | None
) -> Storage | bool:
    for storage in Storage.objects.using(using_key).filter().all():
        # if path is part of storage, return it
        if check_path_is_child_of_root(path, root=storage.root):
            return storage
    return False


def check_path_is_child_of_root(path: Path | UPath, root: Path | UPath | None) -> bool:
    # str is needed to eliminate UPath storage_options
    # from the equality checks below
    path = UPath(str(path))
    root = UPath(str(root))
    return root.resolve() in path.resolve().parents


def get_relative_path_to_directory(
    path: PurePath | Path | UPath, directory: PurePath | Path | UPath
) -> PurePath | Path:
    if isinstance(directory, UPath) and not isinstance(directory, LocalPathClasses):
        # UPath.relative_to() is not behaving as it should (2023-04-07)
        # need to lstrip otherwise inconsistent behavior across trailing slashes
        # see test_artifact.py: test_get_relative_path_to_directory
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


def get_artifact_kwargs_from_data(
    *,
    data: Path | UPath | str | pd.DataFrame | AnnData | MuData,
    key: str | None,
    run: Run | None,
    format: str | None,
    provisional_uid: str,
    default_storage: Storage,
    using_key: str | None = None,
    skip_check_exists: bool = False,
):
    run = get_run(run)
    memory_rep, path, suffix, storage, use_existing_storage_key = process_data(
        provisional_uid,
        data,
        format,
        key,
        default_storage,
        using_key,
        skip_check_exists,
    )
    stat_or_artifact = get_stat_or_artifact(
        path=path,
        using_key=using_key,
    )
    if isinstance(stat_or_artifact, Artifact):
        artifact = stat_or_artifact
        # update the run of the existing artifact
        if run is not None:
            # save the information that this artifact was previously
            # produced by another run
            if artifact.run is not None:
                artifact.run.output_artifacts_with_later_updates.add(artifact)
            # update the run of the artifact with the latest run
            stat_or_artifact.run = run
            stat_or_artifact.transform = run.transform
        return artifact, None
    else:
        size, hash, hash_type, n_objects = stat_or_artifact

    check_path_in_storage = False
    if use_existing_storage_key:
        inferred_key = get_relative_path_to_directory(
            path=path, directory=UPath(storage.root)
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
        storage = default_storage

    if key is not None and key.startswith(AUTO_KEY_PREFIX):
        raise ValueError(f"Key cannot start with {AUTO_KEY_PREFIX}")

    log_storage_hint(
        check_path_in_storage=check_path_in_storage,
        storage=storage,
        key=key,
        uid=provisional_uid,
        suffix=suffix,
        is_dir=n_objects is not None,
    )

    # do we use a virtual or an actual storage key?
    key_is_virtual = settings.creation._artifact_use_virtual_keys

    # if the file is already in storage, independent of the default
    # we use an actual storage key
    if check_path_in_storage:
        key_is_virtual = False

    kwargs = {
        "suffix": suffix,
        "hash": hash,
        "hash_type": hash_type,
        "key": key,
        "size": size,
        "storage_id": storage.id,
        # passing both the id and the object
        # to make them both available immediately
        # after object creation
        "n_objects": n_objects,
        "n_observations": None,  # to implement
        "run_id": run.id if run is not None else None,
        "run": run,
        "key_is_virtual": key_is_virtual,
    }
    if not isinstance(path, LocalPathClasses):
        local_filepath = None
        cloud_filepath = path
    else:
        local_filepath = path
        cloud_filepath = None
    privates = {
        "local_filepath": local_filepath,
        "cloud_filepath": cloud_filepath,
        "memory_rep": memory_rep,
        "check_path_in_storage": check_path_in_storage,
    }
    return kwargs, privates


def log_storage_hint(
    *,
    check_path_in_storage: bool,
    storage: Storage | None,
    key: str | None,
    uid: str,
    suffix: str,
    is_dir: bool,
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
        hint += f"path in storage '{display_root}'"  # type: ignore
    else:
        hint += "path content will be copied to default storage upon `save()`"
    if key is None:
        storage_key = auto_storage_key_from_artifact_uid(uid, suffix, is_dir)
        hint += f" with key `None` ('{storage_key}')"
    else:
        hint += f" with key '{key}'"
    logger.hint(hint)


def data_is_anndata(data: AnnData | UPathStr) -> bool:
    if isinstance(data, AnnData):
        return True
    if isinstance(data, (str, Path, UPath)):
        data_path = UPath(data)
        if data_path.suffix == ".h5ad":
            return True
        elif data_path.suffix == ".zarr":
            # ".anndata.zarr" is a valid suffix (core.storage._valid_suffixes)
            if ".anndata" in data_path.suffixes:
                return True
            # check only for local, expensive for cloud
            if fsspec.utils.get_protocol(data_path.as_posix()) == "file":
                return zarr_is_adata(data_path)
            else:
                logger.warning("We do not check if cloud zarr is AnnData or not.")
                return False
    return False


def data_is_mudata(data: MuData | UPathStr) -> bool:
    if _mudata_is_installed():
        from mudata import MuData

        if isinstance(data, MuData):
            return True
    if isinstance(data, (str, Path)):
        return UPath(data).suffix in {".h5mu"}
    return False


def _check_accessor_artifact(data: Any, accessor: str | None = None):
    if accessor is None:
        if isinstance(data, pd.DataFrame):
            logger.warning("data is a DataFrame, please use .from_df()")
            accessor = "DataFrame"
            return accessor

        data_is_path = isinstance(data, (str, Path))
        if data_is_anndata(data):
            if not data_is_path:
                logger.warning("data is an AnnData, please use .from_anndata()")
            accessor = "AnnData"
        elif data_is_mudata(data):
            if not data_is_path:
                logger.warning("data is a MuData, please use .from_mudata()")
            accessor = "MuData"
        elif not data_is_path:  # UPath is a subclass of Path
            raise TypeError("data has to be a string, Path, UPath")
    return accessor


def update_attributes(data: HasFeatures, attributes: Mapping[str, str]):
    for key, value in attributes.items():
        if getattr(data, key) != value:
            logger.warning(f"updated {key} from {getattr(data, key)} to {value}")
            setattr(data, key, value)


def __init__(artifact: Artifact, *args, **kwargs):
    artifact.features = FeatureManager(artifact)
    artifact.params = ParamManager(artifact)
    # Below checks for the Django-internal call in from_db()
    # it'd be better if we could avoid this, but not being able to create a Artifact
    # from data with the default constructor renders the central class of the API
    # essentially useless
    # The danger below is not that a user might pass as many args (12 of it), but rather
    # that at some point the Django API might change; on the other hand, this
    # condition of for calling the constructor based on kwargs should always
    # stay robust
    if len(args) == len(artifact._meta.concrete_fields):
        super(Artifact, artifact).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 1:
        raise ValueError("Only one non-keyword arg allowed: data")

    data: str | Path = kwargs.pop("data") if len(args) == 0 else args[0]
    type: str = kwargs.pop("type") if "type" in kwargs else "dataset"
    key: str | None = kwargs.pop("key") if "key" in kwargs else None
    run: Run | None = kwargs.pop("run") if "run" in kwargs else None
    description: str | None = (
        kwargs.pop("description") if "description" in kwargs else None
    )
    is_new_version_of: Artifact | None = (
        kwargs.pop("is_new_version_of") if "is_new_version_of" in kwargs else None
    )
    version: str | None = kwargs.pop("version") if "version" in kwargs else None
    visibility: int | None = (
        kwargs.pop("visibility")
        if "visibility" in kwargs
        else VisibilityChoice.default.value
    )
    format = kwargs.pop("format") if "format" in kwargs else None
    skip_check_exists = (
        kwargs.pop("skip_check_exists") if "skip_check_exists" in kwargs else False
    )
    if "default_storage" in kwargs:
        default_storage = kwargs.pop("default_storage")
    else:
        if setup_settings.instance.keep_artifacts_local:
            default_storage = setup_settings.instance.storage_local.record
        else:
            default_storage = setup_settings.instance.storage.record
    using_key = (
        kwargs.pop("using_key") if "using_key" in kwargs else settings._using_key
    )
    accessor = kwargs.pop("accessor") if "accessor" in kwargs else None
    accessor = _check_accessor_artifact(data=data, accessor=accessor)
    if not len(kwargs) == 0:
        raise ValueError(
            "Only data, key, run, description, version, is_new_version_of, visibility"
            f" can be passed, you passed: {kwargs}"
        )

    if is_new_version_of is None:
        provisional_uid = init_uid(version=version, n_full_id=20)
    else:
        if not isinstance(is_new_version_of, Artifact):
            raise TypeError("is_new_version_of has to be of type ln.Artifact")
        provisional_uid, version = get_uid_from_old_version(
            is_new_version_of, version, using_key
        )
        if description is None:
            description = is_new_version_of.description
    kwargs_or_artifact, privates = get_artifact_kwargs_from_data(
        data=data,
        key=key,
        run=run,
        format=format,
        provisional_uid=provisional_uid,
        default_storage=default_storage,
        using_key=using_key,
        skip_check_exists=skip_check_exists,
    )

    # an object with the same hash already exists
    if isinstance(kwargs_or_artifact, Artifact):
        from ._registry import init_self_from_db

        init_self_from_db(artifact, kwargs_or_artifact)
        # adding "key" here is dangerous because key might be auto-populated
        update_attributes(artifact, {"description": description})
        if artifact.key != key and key is not None:
            logger.warning(
                f"key {artifact.key} on existing artifact differs from passed key {key}"
            )
        return None
    else:
        kwargs = kwargs_or_artifact

    # in case we have a new version of a folder with a different hash, print a
    # warning that the old version can't be recovered
    if (
        is_new_version_of is not None
        and is_new_version_of.n_objects is not None
        and is_new_version_of.n_objects > 1
    ):
        logger.warning(
            f"artifact version {version} will _update_ the state of folder {is_new_version_of.path} - "
            "to _retain_ the old state by duplicating the entire folder, do _not_ pass `is_new_version_of`"
        )

    kwargs["type"] = type
    kwargs["uid"] = provisional_uid
    kwargs["version"] = version
    kwargs["description"] = description
    kwargs["visibility"] = visibility
    kwargs["accessor"] = accessor
    # this check needs to come down here because key might be populated from an
    # existing file path during get_artifact_kwargs_from_data()
    if (
        kwargs["key"] is None
        and kwargs["description"] is None
        and kwargs["run"] is None
    ):
        raise ValueError("Pass one of key, run or description as a parameter")

    add_transform_to_kwargs(kwargs, kwargs["run"])

    if data is not None:
        artifact._local_filepath = privates["local_filepath"]
        artifact._cloud_filepath = privates["cloud_filepath"]
        artifact._memory_rep = privates["memory_rep"]
        artifact._to_store = not privates["check_path_in_storage"]

    super(Artifact, artifact).__init__(**kwargs)


@classmethod  # type: ignore
@doc_args(Artifact.from_df.__doc__)
def from_df(
    cls,
    df: pd.DataFrame,
    key: str | None = None,
    description: str | None = None,
    run: Run | None = None,
    version: str | None = None,
    is_new_version_of: Artifact | None = None,
    **kwargs,
) -> Artifact:
    """{}."""
    artifact = Artifact(
        data=df,
        key=key,
        run=run,
        description=description,
        version=version,
        is_new_version_of=is_new_version_of,
        accessor="DataFrame",
        type="dataset",
        **kwargs,
    )
    return artifact


@classmethod  # type: ignore
@doc_args(Artifact.from_anndata.__doc__)
def from_anndata(
    cls,
    adata: AnnData | UPathStr,
    key: str | None = None,
    description: str | None = None,
    run: Run | None = None,
    version: str | None = None,
    is_new_version_of: Artifact | None = None,
    **kwargs,
) -> Artifact:
    """{}."""
    if not data_is_anndata(adata):
        raise ValueError("data has to be an AnnData object or a path to AnnData-like")
    artifact = Artifact(
        data=adata,
        key=key,
        run=run,
        description=description,
        version=version,
        is_new_version_of=is_new_version_of,
        accessor="AnnData",
        type="dataset",
        **kwargs,
    )
    return artifact


@classmethod  # type: ignore
@doc_args(Artifact.from_mudata.__doc__)
def from_mudata(
    cls,
    mdata: MuData,
    key: str | None = None,
    description: str | None = None,
    run: Run | None = None,
    version: str | None = None,
    is_new_version_of: Artifact | None = None,
    **kwargs,
) -> Artifact:
    """{}."""
    artifact = Artifact(
        data=mdata,
        key=key,
        run=run,
        description=description,
        version=version,
        is_new_version_of=is_new_version_of,
        accessor="MuData",
        type="dataset",
        **kwargs,
    )
    return artifact


@classmethod  # type: ignore
@doc_args(Artifact.from_dir.__doc__)
def from_dir(
    cls,
    path: UPathStr,
    key: str | None = None,
    *,
    run: Run | None = None,
) -> list[Artifact]:
    """{}."""
    logger.warning(
        "this creates one artifact per file in the directory - you might simply call"
        " ln.Artifact(dir) to get one artifact for the entire directory"
    )
    folderpath: UPath = create_path(path)  # returns Path for local
    default_storage = settings._storage_settings.record
    using_key = settings._using_key
    storage, use_existing_storage = process_pathlike(
        folderpath, default_storage, using_key
    )
    folder_key_path: PurePath | Path
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
                folderpath, UPath(storage.root)
            )
    else:
        folder_key_path = Path(key)

    # always sanitize by stripping a trailing slash
    folder_key = folder_key_path.as_posix().rstrip("/")

    # TODO: (non-local) UPath doesn't list the first level artifacts and dirs with "*"
    pattern = "" if not isinstance(folderpath, LocalPathClasses) else "*"

    # silence fine-grained logging
    verbosity = settings.verbosity
    verbosity_int = settings._verbosity_int
    if verbosity_int >= 1:
        settings.verbosity = "warning"
    artifacts_dict = {}
    for filepath in folderpath.rglob(pattern):
        if filepath.is_file():
            relative_path = get_relative_path_to_directory(filepath, folderpath)
            artifact_key = folder_key + "/" + relative_path.as_posix()
            # if creating from rglob, we don't need to check for existence
            artifact = Artifact(
                filepath, run=run, key=artifact_key, skip_check_exists=True
            )
            artifacts_dict[artifact.uid] = artifact
    settings.verbosity = verbosity

    # run sanity check on hashes
    hashes = [
        artifact.hash
        for artifact in artifacts_dict.values()
        if artifact.hash is not None
    ]
    uids = artifacts_dict.keys()
    if len(set(hashes)) == len(hashes):
        artifacts = list(artifacts_dict.values())
    else:
        # consider exact duplicates (same id, same hash)
        # below can't happen anymore because artifacts is a dict now
        # if len(set(uids)) == len(set(hashes)):
        #     logger.warning("dropping duplicate records in list of artifact records")
        #     artifacts = list(set(uids))
        # consider false duplicates (different id, same hash)
        if not len(set(uids)) == len(set(hashes)):
            seen_hashes = set()
            non_unique_artifacts = {
                hash: artifact
                for hash, artifact in artifacts_dict.items()
                if artifact.hash in seen_hashes or seen_hashes.add(artifact.hash)  # type: ignore
            }
            display_non_unique = "\n    ".join(
                f"{artifact}" for artifact in non_unique_artifacts
            )
            logger.warning(
                "there are multiple artifact uids with the same hashes, dropping"
                f" {len(non_unique_artifacts)} duplicates out of"
                f" {len(artifacts_dict)} artifacts:\n    {display_non_unique}"
            )
            artifacts = [
                artifact
                for artifact in artifacts_dict.values()
                if artifact not in non_unique_artifacts.values()
            ]
    logger.success(
        f"created {len(artifacts)} artifacts from directory using storage"
        f" {storage.root} and key = {folder_key}/"
    )
    return artifacts


# docstring handled through attach_func_to_class_method
def replace(
    self,
    data: UPathStr,
    run: Run | None = None,
    format: str | None = None,
) -> None:
    default_storage = settings._storage_settings.record
    kwargs, privates = get_artifact_kwargs_from_data(
        provisional_uid=self.uid,
        data=data,
        key=self.key,
        run=run,
        format=format,
        default_storage=default_storage,
    )

    # this artifact already exists
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
        old_storage = auto_storage_key_from_artifact(self)
        is_dir = self.n_objects is not None
        new_storage = auto_storage_key_from_artifact_uid(
            self.uid, kwargs["suffix"], is_dir
        )
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
    self, is_run_input: bool | None = None
) -> AnnDataAccessor | BackedAccessor | SOMACollection | SOMAExperiment:
    # ignore empty suffix for now
    suffixes = (".h5", ".hdf5", ".h5ad", ".zarr", ".tiledbsoma", "")
    if self.suffix not in suffixes:
        raise ValueError(
            "Artifact should have a zarr, h5 or tiledbsoma object as the underlying data, please"
            " use one of the following suffixes for the object name:"
            f" {', '.join(suffixes[:-1])}."
        )

    from lamindb.core.storage._backed_access import backed_access

    _track_run_input(self, is_run_input)
    using_key = settings._using_key
    filepath = filepath_from_artifact(self, using_key=using_key)
    # consider the case where an object is already locally cached
    localpath = setup_settings.instance.storage.cloud_to_local_no_update(filepath)
    if localpath.exists():
        return backed_access(localpath, using_key)
    else:
        return backed_access(filepath, using_key)


# docstring handled through attach_func_to_class_method
def load(self, is_run_input: bool | None = None, stream: bool = False, **kwargs) -> Any:
    _track_run_input(self, is_run_input)
    if hasattr(self, "_memory_rep") and self._memory_rep is not None:
        return self._memory_rep
    using_key = settings._using_key
    return load_to_memory(
        filepath_from_artifact(self, using_key=using_key), stream=stream, **kwargs
    )


# docstring handled through attach_func_to_class_method
def cache(self, is_run_input: bool | None = None) -> Path:
    using_key = settings._using_key
    filepath = filepath_from_artifact(self, using_key=using_key)
    try:
        cache_path = setup_settings.instance.storage.cloud_to_local(
            filepath, print_progress=True
        )
    except Exception as e:
        if not isinstance(filepath, LocalPathClasses):
            cache_path = setup_settings.instance.storage.cloud_to_local_no_update(
                filepath
            )
            if cache_path.is_file():
                cache_path.unlink(missing_ok=True)
            elif cache_path.is_dir():
                shutil.rmtree(cache_path)
        raise e
    # only call if sync is successfull
    _track_run_input(self, is_run_input)
    return cache_path


# docstring handled through attach_func_to_class_method
def delete(
    self,
    permanent: bool | None = None,
    storage: bool | None = None,
    using_key: str | None = None,
) -> None:
    # by default, we only move artifacts into the trash (visibility = -1)
    trash_visibility = VisibilityChoice.trash.value
    if self.visibility > trash_visibility and not permanent:
        if storage is not None:
            logger.warning("moving artifact to trash, storage arg is ignored")
        # move to trash
        self.visibility = trash_visibility
        self.save()
        logger.warning(f"moved artifact to trash (visibility = {trash_visibility})")
        return

    # if the artifact is already in the trash
    # permanent delete skips the trash
    if permanent is None:
        # ask for confirmation of permanent delete
        response = input(
            "Artifact record is already in trash! Are you sure you want to permanently"
            " delete it? (y/n) You can't undo this action."
        )
        delete_record = response == "y"
    else:
        assert permanent  # noqa: S101
        delete_record = True

    if delete_record:
        # need to grab file path before deletion
        try:
            path = filepath_from_artifact(self, using_key)
        except OSError:
            # we can still delete the record
            logger.warning("Could not get path")
            storage = False
        # only delete in storage if DB delete is successful
        # DB delete might error because of a foreign key constraint violated etc.
        self._delete_skip_storage()
        if self.key is None or self.key_is_virtual:
            # do not ask for confirmation also if storage is None
            delete_in_storage = storage is None or storage
        else:
            # for artifacts with non-virtual semantic storage keys (key is not None)
            # ask for extra-confirmation
            if storage is None:
                response = input(
                    f"Are you sure to want to delete {path}? (y/n)  You can't undo"
                    " this action."
                )
                delete_in_storage = response == "y"
            else:
                delete_in_storage = storage
        if not delete_in_storage:
            logger.warning(
                f"you will retain a dangling store here: {path}, not referenced via an artifact"
            )
        # we don't yet have logic to bring back the deleted metadata record
        # in case storage deletion fails - this is important for ACID down the road
        if delete_in_storage:
            delete_msg = delete_storage(path, raise_file_not_found_error=False)
            if delete_msg != "did-not-delete":
                logger.success(f"deleted {colors.yellow(f'{path}')}")


def _delete_skip_storage(artifact, *args, **kwargs) -> None:
    super(Artifact, artifact).delete(*args, **kwargs)


# docstring handled through attach_func_to_class_method
def save(self, upload: bool | None = None, **kwargs) -> None:
    state_was_adding = self._state.adding
    print_progress = kwargs.pop("print_progress", True)
    access_token = kwargs.pop("access_token", None)
    local_path = None
    if upload and setup_settings.instance.keep_artifacts_local:
        # switch local storage location to cloud
        local_path = self.path
        self.storage_id = setup_settings.instance.storage.id
        self._local_filepath = local_path
        # switch to virtual storage key upon upload
        # the local filepath is already cached at that point
        self.key_is_virtual = True
        # ensure that the artifact is uploaded
        self._to_store = True

    self._save_skip_storage(**kwargs)

    from lamindb._save import check_and_attempt_clearing, check_and_attempt_upload

    using_key = None
    if "using" in kwargs:
        using_key = kwargs["using"]
    exception = check_and_attempt_upload(
        self, using_key, access_token=access_token, print_progress=print_progress
    )
    if exception is not None:
        self._delete_skip_storage()
        raise RuntimeError(exception)
    exception = check_and_attempt_clearing(self, using_key)
    if exception is not None:
        raise RuntimeError(exception)
    if local_path is not None and not state_was_adding:
        # only move the local artifact to cache if it was not newly created
        local_path_cache = ln_setup.settings.storage.cache_dir / local_path.name
        # don't use Path.rename here because of cross-device link error
        # https://laminlabs.slack.com/archives/C04A0RMA0SC/p1710259102686969
        shutil.move(
            local_path,  # type: ignore
            local_path_cache,
        )
        logger.important(f"moved local artifact to cache: {local_path_cache}")
    return self


def _save_skip_storage(file, **kwargs) -> None:
    save_feature_sets(file)
    super(Artifact, file).save(**kwargs)
    save_feature_set_links(file)


@property  # type: ignore
@doc_args(Artifact.path.__doc__)
def path(self) -> Path | UPath:
    """{}."""
    using_key = settings._using_key
    return filepath_from_artifact(self, using_key)


# docstring handled through attach_func_to_class_method
def restore(self) -> None:
    self.visibility = VisibilityChoice.default.value
    self.save()


METHOD_NAMES = [
    "__init__",
    "from_anndata",
    "from_df",
    "from_mudata",
    "backed",
    "cache",
    "load",
    "delete",
    "save",
    "replace",
    "from_dir",
    "restore",
]

if ln_setup._TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(Artifact, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, Artifact, globals())

# privates currently dealt with separately
Artifact._delete_skip_storage = _delete_skip_storage
Artifact._save_skip_storage = _save_skip_storage
Artifact.path = path
