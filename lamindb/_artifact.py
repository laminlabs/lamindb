from __future__ import annotations

import os
import shutil
from pathlib import Path, PurePath, PurePosixPath
from typing import TYPE_CHECKING, Any

import fsspec
import lamindb_setup as ln_setup
import pandas as pd
from anndata import AnnData
from django.db.models import Q
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

from lamindb._record import _get_record_kwargs
from lamindb.errors import FieldValidationError
from lamindb.models import Artifact, FeatureManager, ParamManager, Run, Storage

from ._parents import view_lineage
from ._utils import attach_func_to_class_method
from .core._data import (
    _track_run_input,
    describe,
    get_run,
    save_schema_links,
    save_staged_feature_sets,
)
from .core._settings import settings
from .core.loaders import load_to_memory
from .core.storage import (
    LocalPathClasses,
    UPath,
    delete_storage,
    infer_suffix,
    write_to_disk,
)
from .core.storage._anndata_accessor import _anndata_n_observations
from .core.storage._pyarrow_dataset import PYARROW_SUFFIXES
from .core.storage._tiledbsoma import _soma_n_observations
from .core.storage.objects import is_package_installed
from .core.storage.paths import (
    AUTO_KEY_PREFIX,
    auto_storage_key_from_artifact,
    auto_storage_key_from_artifact_uid,
    check_path_is_child_of_root,
    filepath_cache_key_from_artifact,
    filepath_from_artifact,
)
from .core.versioning import (
    create_uid,
    message_update_key_in_version_family,
)
from .errors import IntegrityError, InvalidArgument

try:
    from .core.storage._zarr import identify_zarr_type
except ImportError:

    def identify_zarr_type(storepath):  # type: ignore
        raise ImportError("Please install zarr: pip install zarr<=2.18.4")


if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr
    from mudata import MuData
    from pyarrow.dataset import Dataset as PyArrowDataset
    from spatialdata import SpatialData
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment
    from tiledbsoma import Measurement as SOMAMeasurement

    from lamindb.core.storage._backed_access import AnnDataAccessor, BackedAccessor


def process_pathlike(
    filepath: UPath,
    default_storage: Storage,
    using_key: str | None,
    skip_existence_check: bool = False,
) -> tuple[Storage, bool]:
    """Determines the appropriate storage for a given path and whether to use an existing storage key."""
    if not skip_existence_check:
        try:  # check if file exists
            if not filepath.exists():
                raise FileNotFoundError(filepath)
        except PermissionError:
            pass
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
                if filepath.protocol == "hf":
                    hf_path = filepath.fs.resolve_path(filepath.as_posix())
                    hf_path.path_in_repo = ""
                    new_root = "hf://" + hf_path.unresolve()
                else:
                    if filepath.protocol == "s3":
                        # check that endpoint_url didn't propagate here
                        # as a part of the path string
                        assert "?" not in filepath.path  # noqa: S101
                    new_root = list(filepath.parents)[-1]
                # do not register remote storage locations on hub if the current instance
                # is not managed on the hub
                storage_settings, _ = init_storage(
                    new_root, prevent_register_hub=not setup_settings.instance.is_on_hub
                )
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
    is_replace: bool = False,
) -> tuple[Any, Path | UPath, str, Storage, bool]:
    """Serialize a data object that's provided as file or in memory.

    if not overwritten, data gets stored in default storage
    """
    supported_data_types = [pd.DataFrame, AnnData]
    if is_package_installed("mudata"):
        from mudata import MuData

        supported_data_types.append(MuData)
    if is_package_installed("spatialdata"):
        from spatialdata import SpatialData

        supported_data_types.append(SpatialData)
    supported_data_types = tuple(supported_data_types)  # type: ignore

    if key is not None:
        key_suffix = extract_suffix_from_path(PurePosixPath(key), arg_name="key")
        # use suffix as the (adata) format if the format is not provided
        if isinstance(data, AnnData) and format is None and len(key_suffix) > 0:
            format = key_suffix[1:]
    else:
        key_suffix = None
    if isinstance(data, (str, Path, UPath)):  # UPathStr, spelled out
        access_token = (
            default_storage._access_token
            if hasattr(default_storage, "_access_token")
            else None
        )
        path = create_path(data, access_token=access_token)
        # we don't resolve http links because they can resolve into a different domain
        # for example into a temporary url
        if path.protocol not in {"http", "https"}:
            path = path.resolve()
        storage, use_existing_storage_key = process_pathlike(
            path,
            default_storage=default_storage,
            using_key=using_key,
            skip_existence_check=skip_existence_check,
        )
        suffix = extract_suffix_from_path(path)
        memory_rep = None
    elif isinstance(data, supported_data_types):
        storage = default_storage
        memory_rep = data
        suffix = infer_suffix(data, format)
    else:
        raise NotImplementedError(
            f"Do not know how to create a artifact object from {data}, pass a path instead!"
        )
    if key_suffix is not None and key_suffix != suffix and not is_replace:
        # consciously omitting a trailing period
        if isinstance(data, (str, Path, UPath)):
            message = f"The suffix '{suffix}' of the provided path is inconsistent, it should be '{key_suffix}'"
        else:
            message = f"The suffix '{key_suffix}' of the provided key is inconsistent, it should be '{suffix}'"
        raise InvalidArgument(message)
    # in case we have an in-memory representation, we need to write it to disk
    if isinstance(data, supported_data_types):
        path = settings.cache_dir / f"{provisional_uid}{suffix}"
        write_to_disk(data, path)
        use_existing_storage_key = False
    return memory_rep, path, suffix, storage, use_existing_storage_key


def get_stat_or_artifact(
    path: UPath,
    key: str | None = None,
    check_hash: bool = True,
    is_replace: bool = False,
    instance: str | None = None,
) -> tuple[int, str | None, str | None, int | None, Artifact | None] | Artifact:
    """Retrieves file statistics or an existing artifact based on the path, hash, and key."""
    n_files = None
    if settings.creation.artifact_skip_size_hash:
        return None, None, None, n_files, None
    stat = path.stat()  # one network request
    if not isinstance(path, LocalPathClasses):
        size, hash, hash_type = None, None, None
        if stat is not None:
            # convert UPathStatResult to fsspec info dict
            stat = stat.as_info()
            if (store_type := stat["type"]) == "file":
                size, hash, hash_type = get_stat_file_cloud(stat)
            elif store_type == "directory":
                size, hash, hash_type, n_files = get_stat_dir_cloud(path)
        if hash is None:
            logger.warning(f"did not add hash for {path}")
            return size, hash, hash_type, n_files, None
    else:
        if path.is_dir():
            size, hash, hash_type, n_files = hash_dir(path)
        else:
            hash, hash_type = hash_file(path)
            size = stat.st_size
    if not check_hash:
        return size, hash, hash_type, n_files, None
    previous_artifact_version = None
    if key is None or is_replace:
        result = Artifact.objects.using(instance).filter(hash=hash).all()
        artifact_with_same_hash_exists = len(result) > 0
    else:
        storage_id = settings.storage.id
        result = (
            Artifact.objects.using(instance)
            .filter(Q(hash=hash) | Q(key=key, storage_id=storage_id))
            .order_by("-created_at")
            .all()
        )
        artifact_with_same_hash_exists = result.filter(hash=hash).count() > 0
        if not artifact_with_same_hash_exists and len(result) > 0:
            logger.important(
                f"creating new artifact version for key='{key}' (storage: '{settings.storage.root_as_str}')"
            )
            previous_artifact_version = result[0]
    if artifact_with_same_hash_exists:
        message = "found artifact with same hash"
        if result[0]._branch_code == -1:
            result[0].restore()
            message = "restored artifact with same hash from trash"
        logger.important(
            f"{message}: {result[0]}; to track this artifact as an input, use: ln.Artifact.get()"
        )
        return result[0]
    else:
        return size, hash, hash_type, n_files, previous_artifact_version


def check_path_in_existing_storage(
    path: Path | UPath, using_key: str | None = None
) -> Storage | bool:
    for storage in Storage.objects.using(using_key).filter().all():
        # if path is part of storage, return it
        if check_path_is_child_of_root(path, root=storage.root):
            return storage
    return False


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
    version: str | None,
    default_storage: Storage,
    using_key: str | None = None,
    is_replace: bool = False,
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
        is_replace=is_replace,
    )
    stat_or_artifact = get_stat_or_artifact(
        path=path,
        key=key,
        instance=using_key,
        is_replace=is_replace,
    )
    if isinstance(stat_or_artifact, Artifact):
        artifact = stat_or_artifact
        # update the run of the existing artifact
        if run is not None:
            # save the information that this artifact was previously produced by
            # another run
            # note: same logic exists for _output_collections_with_later_updates
            if artifact.run is not None and artifact.run != run:
                artifact.run._output_artifacts_with_later_updates.add(artifact)
            # update the run of the artifact with the latest run
            stat_or_artifact.run = run
        return artifact, None
    else:
        size, hash, hash_type, n_files, revises = stat_or_artifact

    if revises is not None:  # update provisional_uid
        provisional_uid, revises = create_uid(revises=revises, version=version)
        if settings.cache_dir in path.parents:
            path = path.rename(path.with_name(f"{provisional_uid}{suffix}"))

    check_path_in_storage = False
    if use_existing_storage_key:
        inferred_key = get_relative_path_to_directory(
            path=path, directory=UPath(storage.root)
        ).as_posix()
        if key is None:
            key = inferred_key
        else:
            if not key == inferred_key:
                raise InvalidArgument(
                    f"The path '{data}' is already in registered storage"
                    f" '{storage.root}' with key '{inferred_key}'\nYou passed"
                    f" conflicting key '{key}': please move the file before"
                    " registering it."
                )
        check_path_in_storage = True
    else:
        storage = default_storage

    log_storage_hint(
        check_path_in_storage=check_path_in_storage,
        storage=storage,
        key=key,
        uid=provisional_uid,
        suffix=suffix,
        is_dir=n_files is not None,
    )

    # do we use a virtual or an actual storage key?
    key_is_virtual = settings.creation._artifact_use_virtual_keys

    # if the file is already in storage, independent of the default
    # we use an actual storage key
    if check_path_in_storage:
        key_is_virtual = False

    kwargs = {
        "uid": provisional_uid,
        "suffix": suffix,
        "hash": hash,
        "_hash_type": hash_type,
        "key": key,
        "size": size,
        "storage_id": storage.id,
        # passing both the id and the object
        # to make them both available immediately
        # after object creation
        "n_files": n_files,
        "_overwrite_versions": n_files is not None,  # True for folder, False for file
        "n_observations": None,  # to implement
        "run_id": run.id if run is not None else None,
        "run": run,
        "_key_is_virtual": key_is_virtual,
        "revises": revises,
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
                display_root = root_path.relative_to(Path.cwd())  # type: ignore
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
        if ".h5ad" in data_path.suffixes:  # ".h5ad.gz" is a valid suffix
            return True
        elif data_path.suffix == ".zarr":
            # ".anndata.zarr" is a valid suffix (core.storage._valid_suffixes)
            if ".anndata" in data_path.suffixes:
                return True
            # check only for local, expensive for cloud
            if fsspec.utils.get_protocol(data_path.as_posix()) == "file":
                return identify_zarr_type(data_path) == "anndata"
            else:
                logger.warning("We do not check if cloud zarr is AnnData or not")
                return False
    return False


def data_is_mudata(data: MuData | UPathStr) -> bool:
    if is_package_installed("mudata"):
        from mudata import MuData

        if isinstance(data, MuData):
            return True
    if isinstance(data, (str, Path)):
        return UPath(data).suffix == ".h5mu"
    return False


def data_is_spatialdata(data: SpatialData | UPathStr) -> bool:
    if is_package_installed("spatialdata"):
        from spatialdata import SpatialData

        if isinstance(data, SpatialData):
            return True
        if isinstance(data, (str, Path)):
            if UPath(data).suffix in ["spatialdata.zarr", ".zarr"]:
                return identify_zarr_type(data, check=False) == "spatialdata"
        return False


def _check_otype_artifact(
    data: UPathStr | pd.DataFrame | AnnData | MuData | SpatialData,
    otype: str | None = None,
) -> str:
    if otype is None:
        if isinstance(data, pd.DataFrame):
            logger.warning("data is a DataFrame, please use .from_df()")
            otype = "DataFrame"
            return otype

        data_is_path = isinstance(data, (str, Path))
        if data_is_anndata(data):
            if not data_is_path:
                logger.warning("data is an AnnData, please use .from_anndata()")
            otype = "AnnData"
        elif data_is_mudata(data):
            if not data_is_path:
                logger.warning("data is a MuData, please use .from_mudata()")
            otype = "MuData"
        elif data_is_spatialdata(data):
            if not data_is_path:
                logger.warning("data is a SpatialData, please use .from_spatialdata()")
            otype = "SpatialData"
        elif not data_is_path:  # UPath is a subclass of Path
            raise TypeError("data has to be a string, Path, UPath")
    return otype


def __init__(artifact: Artifact, *args, **kwargs):
    artifact.features = FeatureManager(artifact)  # type: ignore
    artifact.params = ParamManager(artifact)  # type: ignore
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
    kind: str = kwargs.pop("kind") if "kind" in kwargs else None
    key: str | None = kwargs.pop("key") if "key" in kwargs else None
    run: Run | None = kwargs.pop("run") if "run" in kwargs else None
    description: str | None = (
        kwargs.pop("description") if "description" in kwargs else None
    )
    revises: Artifact | None = kwargs.pop("revises") if "revises" in kwargs else None
    version: str | None = kwargs.pop("version") if "version" in kwargs else None
    if "visibility" in kwargs:
        _branch_code = kwargs.pop("visibility")
    elif "_branch_code" in kwargs:
        _branch_code = kwargs.pop("_branch_code")
    else:
        _branch_code = 1
    format = kwargs.pop("format") if "format" in kwargs else None
    _is_internal_call = kwargs.pop("_is_internal_call", False)
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
    otype = kwargs.pop("otype") if "otype" in kwargs else None
    otype = _check_otype_artifact(data=data, otype=otype)
    if "type" in kwargs:
        logger.warning("`type` will be removed soon, please use `kind`")
        kind = kwargs.pop("type")
    if not len(kwargs) == 0:
        valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Artifact)])
        raise FieldValidationError(
            f"Only {valid_keywords} can be passed, you passed: {kwargs}"
        )
    if revises is not None and key is not None and revises.key != key:
        note = message_update_key_in_version_family(
            suid=revises.stem_uid,
            existing_key=revises.key,
            new_key=key,
            registry="Artifact",
        )
        raise ValueError(
            f"`key` is {key}, but `revises.key` is '{revises.key}'\n\n Either do *not* pass `key`.\n\n{note}"
        )
    if revises is not None:
        if not isinstance(revises, Artifact):
            raise TypeError("`revises` has to be of type `Artifact`")
        if description is None:
            description = revises.description
    if key is not None and AUTO_KEY_PREFIX in key:
        raise ValueError(
            f"Do not pass key that contains a managed storage path in `{AUTO_KEY_PREFIX}`"
        )
    # below is for internal calls that require defining the storage location
    # ahead of constructing the Artifact
    if isinstance(data, (str, Path)) and AUTO_KEY_PREFIX in str(data):
        if _is_internal_call:
            is_automanaged_path = True
            user_provided_key = key
            key = None
        else:
            raise ValueError(
                f"Do not pass path inside the `{AUTO_KEY_PREFIX}` directory."
            )
    else:
        is_automanaged_path = False
    provisional_uid, revises = create_uid(revises=revises, version=version)
    kwargs_or_artifact, privates = get_artifact_kwargs_from_data(
        data=data,
        key=key,
        run=run,
        format=format,
        provisional_uid=provisional_uid,
        version=version,
        default_storage=default_storage,
        using_key=using_key,
        skip_check_exists=skip_check_exists,
    )

    # an object with the same hash already exists
    if isinstance(kwargs_or_artifact, Artifact):
        from ._record import init_self_from_db, update_attributes

        init_self_from_db(artifact, kwargs_or_artifact)
        # adding "key" here is dangerous because key might be auto-populated
        attr_to_update = {"description": description}
        if kwargs_or_artifact._key_is_virtual and kwargs_or_artifact.key is None:
            attr_to_update["key"] = key
        elif artifact.key != key and key is not None:
            logger.warning(
                f"key {artifact.key} on existing artifact differs from passed key {key}"
            )
        update_attributes(artifact, attr_to_update)
        return None
    else:
        kwargs = kwargs_or_artifact

    if revises is None:
        revises = kwargs_or_artifact.pop("revises")

    if data is not None:
        artifact._local_filepath = privates["local_filepath"]
        artifact._cloud_filepath = privates["cloud_filepath"]
        artifact._memory_rep = privates["memory_rep"]
        artifact._to_store = not privates["check_path_in_storage"]

    if is_automanaged_path and _is_internal_call:
        kwargs["_key_is_virtual"] = True
        assert AUTO_KEY_PREFIX in kwargs["key"]  # noqa: S101
        uid = kwargs["key"].replace(AUTO_KEY_PREFIX, "").replace(kwargs["suffix"], "")
        kwargs["key"] = user_provided_key
        if revises is not None:
            assert uid.startswith(revises.stem_uid)  # noqa: S101
        if len(uid) == 16:
            if revises is None:
                uid += "0000"
            else:
                uid, revises = create_uid(revises=revises, version=version)
        kwargs["uid"] = uid

    # only set key now so that we don't do a look-up on it in case revises is passed
    if revises is not None:
        kwargs["key"] = revises.key

    kwargs["kind"] = kind
    kwargs["version"] = version
    kwargs["description"] = description
    kwargs["_branch_code"] = _branch_code
    kwargs["otype"] = otype
    kwargs["revises"] = revises
    # this check needs to come down here because key might be populated from an
    # existing file path during get_artifact_kwargs_from_data()
    if (
        kwargs["key"] is None
        and kwargs["description"] is None
        and kwargs["run"] is None
    ):
        raise ValueError("Pass one of key, run or description as a parameter")

    super(Artifact, artifact).__init__(**kwargs)


@classmethod  # type: ignore
@doc_args(Artifact.from_df.__doc__)
def from_df(
    cls,
    df: pd.DataFrame,
    *,
    key: str | None = None,
    description: str | None = None,
    run: Run | None = None,
    revises: Artifact | None = None,
    **kwargs,
) -> Artifact:
    """{}"""  # noqa: D415
    artifact = Artifact(  # type: ignore
        data=df,
        key=key,
        run=run,
        description=description,
        revises=revises,
        otype="DataFrame",
        kind="dataset",
        **kwargs,
    )
    artifact.n_observations = len(df)
    return artifact


@classmethod  # type: ignore
@doc_args(Artifact.from_anndata.__doc__)
def from_anndata(
    cls,
    adata: AnnData | UPathStr,
    *,
    key: str | None = None,
    description: str | None = None,
    run: Run | None = None,
    revises: Artifact | None = None,
    **kwargs,
) -> Artifact:
    """{}"""  # noqa: D415
    if not data_is_anndata(adata):
        raise ValueError("data has to be an AnnData object or a path to AnnData-like")
    _anndata_n_observations(adata)
    artifact = Artifact(  # type: ignore
        data=adata,
        key=key,
        run=run,
        description=description,
        revises=revises,
        otype="AnnData",
        kind="dataset",
        **kwargs,
    )
    # this is done instead of _anndata_n_observations(adata)
    # because we need a proper path through create_path for cloud paths
    # for additional upath options etc that create_path adds
    obj_for_obs: AnnData | UPath
    if hasattr(artifact, "_memory_rep") and artifact._memory_rep is not None:
        obj_for_obs = artifact._memory_rep
    else:
        # returns ._local_filepath for local files
        # and the proper path through create_path for cloud paths
        obj_for_obs = artifact.path
    artifact.n_observations = _anndata_n_observations(obj_for_obs)
    return artifact


@classmethod  # type: ignore
@doc_args(Artifact.from_mudata.__doc__)
def from_mudata(
    cls,
    mdata: MuData | UPathStr,
    *,
    key: str | None = None,
    description: str | None = None,
    run: Run | None = None,
    revises: Artifact | None = None,
    **kwargs,
) -> Artifact:
    """{}"""  # noqa: D415
    if not data_is_mudata(mdata):
        raise ValueError("data has to be an MuData object or a path to MuData-like")
    artifact = Artifact(  # type: ignore
        data=mdata,
        key=key,
        run=run,
        description=description,
        revises=revises,
        otype="MuData",
        kind="dataset",
        **kwargs,
    )
    artifact.n_observations = mdata.n_obs
    return artifact


@classmethod  # type: ignore
@doc_args(Artifact.from_spatialdata.__doc__)
def from_spatialdata(
    cls,
    sdata: SpatialData | UPathStr,
    *,
    key: str | None = None,
    description: str | None = None,
    run: Run | None = None,
    revises: Artifact | None = None,
    **kwargs,
) -> Artifact:
    """{}"""  # noqa: D415
    if not data_is_spatialdata(sdata):
        raise ValueError("data has to be an MuData object or a path to MuData-like")
    artifact = Artifact(  # type: ignore
        data=sdata,
        key=key,
        run=run,
        description=description,
        revises=revises,
        otype="SpatialData",
        kind="dataset",
        **kwargs,
    )
    # ill-defined https://scverse.zulipchat.com/#narrow/channel/315824-spatial/topic/How.20to.20calculate.20the.20number.20of.20observations.3F
    # artifact.n_observations = ...
    return artifact


@classmethod  # type: ignore
@doc_args(Artifact.from_tiledbsoma.__doc__)
def from_tiledbsoma(
    cls,
    path: UPathStr,
    *,
    key: str | None = None,
    description: str | None = None,
    run: Run | None = None,
    revises: Artifact | None = None,
    **kwargs,
) -> Artifact:
    """{}"""  # noqa: D415
    if UPath(path).suffix != ".tiledbsoma":
        raise ValueError(
            "A tiledbsoma store should have .tiledbsoma suffix to be registered."
        )
    artifact = Artifact(  # type: ignore
        data=path,
        key=key,
        run=run,
        description=description,
        revises=revises,
        otype="tiledbsoma",
        kind="dataset",
        **kwargs,
    )
    artifact.n_observations = _soma_n_observations(artifact.path)
    return artifact


@classmethod  # type: ignore
@doc_args(Artifact.from_dir.__doc__)
def from_dir(
    cls,
    path: UPathStr,
    *,
    key: str | None = None,
    run: Run | None = None,
) -> list[Artifact]:
    """{}"""  # noqa: D415
    folderpath: UPath = create_path(path)  # returns Path for local
    default_storage = settings.storage.record
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

    folder_key = folder_key_path.as_posix()
    # silence fine-grained logging
    verbosity = settings.verbosity
    verbosity_int = settings._verbosity_int
    if verbosity_int >= 1:
        settings.verbosity = "warning"
    artifacts_dict = {}
    for filepath in folderpath.rglob("*"):
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
    n_unique_hashes = len(set(hashes))
    if n_unique_hashes == len(hashes):
        artifacts = list(artifacts_dict.values())
    else:
        # consider exact duplicates (same id, same hash)
        # below can't happen anymore because artifacts is a dict now
        # if len(set(uids)) == len(set(hashes)):
        #     logger.warning("dropping duplicate records in list of artifact records")
        #     artifacts = list(set(uids))
        # consider false duplicates (different id, same hash)
        if not len(set(uids)) == n_unique_hashes:
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
    data: UPathStr | pd.DataFrame | AnnData | MuData,
    run: Run | None = None,
    format: str | None = None,
) -> None:
    default_storage = settings.storage.record
    kwargs, privates = get_artifact_kwargs_from_data(
        provisional_uid=self.uid,
        data=data,
        key=self.key,
        run=run,
        format=format,
        default_storage=default_storage,
        version=None,
        is_replace=True,
    )

    # this artifact already exists
    if privates is None:
        return kwargs

    check_path_in_storage = privates["check_path_in_storage"]
    if check_path_in_storage:
        err_msg = (
            "Can only replace with a local path not in any Storage. "
            f"This data is in {Storage.objects.get(id=kwargs['storage_id'])}."
        )
        raise ValueError(err_msg)

    _overwrite_versions = kwargs["_overwrite_versions"]
    if self._overwrite_versions != _overwrite_versions:
        err_msg = "It is not allowed to replace "
        err_msg += "a folder" if self._overwrite_versions else "a file"
        err_msg += " with " + ("a folder." if _overwrite_versions else "a file.")
        raise ValueError(err_msg)

    if self.key is not None and not self._key_is_virtual:
        key_path = PurePosixPath(self.key)
        new_filename = f"{key_path.stem}{kwargs['suffix']}"
        # the following will only be true if the suffix changes!
        if key_path.name != new_filename:
            self._clear_storagekey = self.key
            self.key = str(key_path.with_name(new_filename))
            # update old key with the new one so that checks in record pass
            self._old_key = self.key
            logger.warning(
                f"replacing the file will replace key '{key_path}' with '{self.key}'"
                f" and delete '{key_path}' upon `save()`"
            )
    else:
        old_storage = auto_storage_key_from_artifact(self)
        is_dir = self.n_files is not None
        new_storage = auto_storage_key_from_artifact_uid(
            self.uid, kwargs["suffix"], is_dir
        )
        if old_storage != new_storage:
            self._clear_storagekey = old_storage
            if self.key is not None:
                new_key_path = PurePosixPath(self.key).with_suffix(kwargs["suffix"])
                self.key = str(new_key_path)
                # update old key with the new one so that checks in record pass
                self._old_key = self.key

    self.suffix = kwargs["suffix"]
    self.size = kwargs["size"]
    self.hash = kwargs["hash"]
    self._hash_type = kwargs["_hash_type"]
    self.run_id = kwargs["run_id"]
    self.run = kwargs["run"]
    self.n_files = kwargs["n_files"]

    self._local_filepath = privates["local_filepath"]
    self._cloud_filepath = privates["cloud_filepath"]
    self._memory_rep = privates["memory_rep"]
    # no need to upload if new file is already in storage
    self._to_store = not check_path_in_storage


inconsistent_state_msg = (
    "Trying to read a folder artifact from an outdated version, "
    "this can result in an incosistent state.\n"
    "Read from the latest version: artifact.versions.filter(is_latest=True).one()"
)


# docstring handled through attach_func_to_class_method
def open(
    self, mode: str = "r", is_run_input: bool | None = None, **kwargs
) -> (
    AnnDataAccessor
    | BackedAccessor
    | SOMACollection
    | SOMAExperiment
    | SOMAMeasurement
    | PyArrowDataset
):
    if self._overwrite_versions and not self.is_latest:
        raise ValueError(inconsistent_state_msg)
    # all hdf5 suffixes including gzipped
    h5_suffixes = [".h5", ".hdf5", ".h5ad"]
    h5_suffixes += [s + ".gz" for s in h5_suffixes]
    # ignore empty suffix for now
    suffixes = (
        (
            "",
            ".zarr",
            ".anndata.zarr",
            ".tiledbsoma",
        )
        + tuple(h5_suffixes)
        + PYARROW_SUFFIXES
        + tuple(
            s + ".gz" for s in PYARROW_SUFFIXES
        )  # this doesn't work for externally gzipped files, REMOVE LATER
    )
    if self.suffix not in suffixes:
        raise ValueError(
            "Artifact should have a zarr, h5, tiledbsoma object"
            " or a compatible `pyarrow.dataset.dataset` directory"
            " as the underlying data, please use one of the following suffixes"
            f" for the object name: {', '.join(suffixes[1:])}."
            f" Or no suffix for a folder with {', '.join(PYARROW_SUFFIXES)} files"
            " (no mixing allowed)."
        )
    if self.suffix != ".tiledbsoma" and self.key != "soma" and mode != "r":
        raise ValueError("Only a tiledbsoma store can be openened with `mode!='r'`.")

    from lamindb.core.storage._backed_access import _track_writes_factory, backed_access

    using_key = settings._using_key
    filepath, cache_key = filepath_cache_key_from_artifact(self, using_key=using_key)
    is_tiledbsoma_w = (
        filepath.name == "soma" or self.suffix == ".tiledbsoma"
    ) and mode == "w"
    # consider the case where an object is already locally cached
    localpath = setup_settings.paths.cloud_to_local_no_update(
        filepath, cache_key=cache_key
    )
    if is_tiledbsoma_w:
        open_cache = False
    else:
        open_cache = not isinstance(
            filepath, LocalPathClasses
        ) and not filepath.synchronize(localpath, just_check=True)
    if open_cache:
        try:
            access = backed_access(localpath, mode, using_key, **kwargs)
        except Exception as e:
            if isinstance(filepath, LocalPathClasses):
                raise e
            logger.warning(
                f"The cache might be corrupted: {e}. Trying to open directly."
            )
            access = backed_access(filepath, mode, using_key, **kwargs)
            # happens only if backed_access has been successful
            # delete the corrupted cache
            if localpath.is_dir():
                shutil.rmtree(localpath)
            else:
                localpath.unlink(missing_ok=True)
    else:
        access = backed_access(filepath, mode, using_key, **kwargs)
        if is_tiledbsoma_w:

            def finalize():
                nonlocal self, filepath, localpath
                if not isinstance(filepath, LocalPathClasses):
                    _, hash, _, _ = get_stat_dir_cloud(filepath)
                else:
                    # this can be very slow
                    _, hash, _, _ = hash_dir(filepath)
                if self.hash != hash:
                    from ._record import init_self_from_db

                    new_version = Artifact(
                        filepath, revises=self, _is_internal_call=True
                    ).save()
                    init_self_from_db(self, new_version)

                    if localpath != filepath and localpath.exists():
                        shutil.rmtree(localpath)

            access = _track_writes_factory(access, finalize)
    # only call if open is successfull
    _track_run_input(self, is_run_input)
    return access


# can't really just call .cache in .load because of double tracking
def _synchronize_cleanup_on_error(
    filepath: UPath, cache_key: str | None = None
) -> UPath:
    try:
        cache_path = setup_settings.paths.cloud_to_local(
            filepath, cache_key=cache_key, print_progress=True
        )
    except Exception as e:
        if not isinstance(filepath, LocalPathClasses):
            cache_path = setup_settings.paths.cloud_to_local_no_update(
                filepath, cache_key=cache_key
            )
            if cache_path.is_dir():
                shutil.rmtree(cache_path)
            else:
                cache_path.unlink(missing_ok=True)
        raise e
    return cache_path


# docstring handled through attach_func_to_class_method
def load(self, is_run_input: bool | None = None, **kwargs) -> Any:
    if self._overwrite_versions and not self.is_latest:
        raise ValueError(inconsistent_state_msg)

    if hasattr(self, "_memory_rep") and self._memory_rep is not None:
        access_memory = self._memory_rep
    else:
        filepath, cache_key = filepath_cache_key_from_artifact(
            self, using_key=settings._using_key
        )
        cache_path = _synchronize_cleanup_on_error(filepath, cache_key=cache_key)
        try:
            # cache_path is local so doesn't trigger any sync in load_to_memory
            access_memory = load_to_memory(cache_path, **kwargs)
        except Exception as e:
            # just raise the exception if the original path is local
            if isinstance(filepath, LocalPathClasses):
                raise e
            logger.warning(
                f"The cache might be corrupted: {e}. Retrying to synchronize."
            )
            # delete the existing cache
            if cache_path.is_dir():
                shutil.rmtree(cache_path)
            else:
                cache_path.unlink(missing_ok=True)
            # download again and try to load into memory
            cache_path = _synchronize_cleanup_on_error(filepath, cache_key=cache_key)
            access_memory = load_to_memory(cache_path, **kwargs)
    # only call if load is successfull
    _track_run_input(self, is_run_input)
    return access_memory


# docstring handled through attach_func_to_class_method
def cache(self, is_run_input: bool | None = None) -> Path:
    if self._overwrite_versions and not self.is_latest:
        raise ValueError(inconsistent_state_msg)

    filepath, cache_key = filepath_cache_key_from_artifact(
        self, using_key=settings._using_key
    )
    cache_path = _synchronize_cleanup_on_error(filepath, cache_key=cache_key)
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
    # this first check means an invalid delete fails fast rather than cascading through
    # database and storage permission errors
    if os.getenv("LAMINDB_MULTI_INSTANCE") is None:
        isettings = setup_settings.instance
        if self.storage.instance_uid != isettings.uid and (storage or storage is None):
            raise IntegrityError(
                "Cannot simply delete artifacts outside of this instance's managed storage locations."
                "\n(1) If you only want to delete the metadata record in this instance, pass `storage=False`"
                f"\n(2) If you want to delete the artifact in storage, please load the managing lamindb instance (uid={self.storage.instance_uid})."
                f"\nThese are all managed storage locations of this instance:\n{Storage.filter(instance_uid=isettings.uid).df()}"
            )
    # by default, we only move artifacts into the trash (_branch_code = -1)
    trash__branch_code = -1
    if self._branch_code > trash__branch_code and not permanent:
        if storage is not None:
            logger.warning("moving artifact to trash, storage arg is ignored")
        # move to trash
        self._branch_code = trash__branch_code
        self.save()
        logger.important(
            f"moved artifact to trash (_branch_code = {trash__branch_code})"
        )
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
            path, _ = filepath_from_artifact(self, using_key)
        except OSError:
            # we can still delete the record
            logger.warning("Could not get path")
            storage = False
        # only delete in storage if DB delete is successful
        # DB delete might error because of a foreign key constraint violated etc.
        if self._overwrite_versions and self.is_latest:
            # includes self
            for version in self.versions.all():
                _delete_skip_storage(version)
        else:
            self._delete_skip_storage()
        # by default do not delete storage if deleting only a previous version
        # and the underlying store is mutable
        if self._overwrite_versions and not self.is_latest:
            delete_in_storage = False
            if storage:
                logger.warning(
                    "Storage argument is ignored; can't delete storage on an previous version"
                )
        elif self.key is None or self._key_is_virtual:
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
            logger.important(f"a file/folder remains here: {path}")
        # we don't yet have logic to bring back the deleted metadata record
        # in case storage deletion fails - this is important for ACID down the road
        if delete_in_storage:
            delete_msg = delete_storage(path, raise_file_not_found_error=False)
            if delete_msg != "did-not-delete":
                logger.success(f"deleted {colors.yellow(f'{path}')}")


def _delete_skip_storage(artifact, *args, **kwargs) -> None:
    super(Artifact, artifact).delete(*args, **kwargs)


# docstring handled through attach_func_to_class_method
def save(self, upload: bool | None = None, **kwargs) -> Artifact:
    state_was_adding = self._state.adding
    print_progress = kwargs.pop("print_progress", True)
    store_kwargs = kwargs.pop("store_kwargs", {})  # kwargs for .upload_from in the end
    access_token = kwargs.pop("access_token", None)
    local_path = None
    if upload and setup_settings.instance.keep_artifacts_local:
        # switch local storage location to cloud
        local_path = self.path
        self.storage_id = setup_settings.instance.storage.id
        self._local_filepath = local_path
        # switch to virtual storage key upon upload
        # the local filepath is already cached at that point
        self._key_is_virtual = True
        # ensure that the artifact is uploaded
        self._to_store = True

    self._save_skip_storage(**kwargs)

    from lamindb._save import check_and_attempt_clearing, check_and_attempt_upload

    using_key = None
    if "using" in kwargs:
        using_key = kwargs["using"]
    exception_upload = check_and_attempt_upload(
        self,
        using_key,
        access_token=access_token,
        print_progress=print_progress,
        **store_kwargs,
    )
    if exception_upload is not None:
        # we do not want to raise file not found on cleanup if upload of a file failed
        # often it is ACID in the filesystem itself
        # for example, s3 won't have the failed file, so just skip the delete in this case
        raise_file_not_found_error = False
        self._delete_skip_storage()
    else:
        # this is the case when it is cleaned on .replace
        raise_file_not_found_error = True
    # this is triggered by an exception in check_and_attempt_upload or by replace.
    exception_clear = check_and_attempt_clearing(
        self, raise_file_not_found_error=raise_file_not_found_error, using_key=using_key
    )
    if exception_upload is not None:
        raise RuntimeError(exception_upload)
    if exception_clear is not None:
        raise RuntimeError(exception_clear)
    # this is only for keep_artifacts_local
    if local_path is not None and not state_was_adding:
        # only move the local artifact to cache if it was not newly created
        local_path_cache = ln_setup.settings.cache_dir / local_path.name
        # don't use Path.rename here because of cross-device link error
        # https://laminlabs.slack.com/archives/C04A0RMA0SC/p1710259102686969
        shutil.move(
            local_path,  # type: ignore
            local_path_cache,
        )
        logger.important(f"moved local artifact to cache: {local_path_cache}")
    return self


def _save_skip_storage(file, **kwargs) -> None:
    save_staged_feature_sets(file)
    super(Artifact, file).save(**kwargs)
    save_schema_links(file)


@property  # type: ignore
@doc_args(Artifact.path.__doc__)
def path(self) -> Path | UPath:
    """{}"""  # noqa: D415
    # return only the path, without StorageSettings
    filepath, _ = filepath_from_artifact(self, using_key=settings._using_key)
    return filepath


# get cache path without triggering sync
@property  # type: ignore
def _cache_path(self) -> UPath:
    filepath, cache_key = filepath_cache_key_from_artifact(
        self, using_key=settings._using_key
    )
    if isinstance(filepath, LocalPathClasses):
        return filepath
    return setup_settings.paths.cloud_to_local_no_update(filepath, cache_key=cache_key)


# docstring handled through attach_func_to_class_method
def restore(self) -> None:
    self._branch_code = 1
    self.save()


METHOD_NAMES = [
    "__init__",
    "from_anndata",
    "from_df",
    "from_mudata",
    "from_spatialdata",
    "from_tiledbsoma",
    "open",
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
# mypy: ignore-errors
Artifact._delete_skip_storage = _delete_skip_storage
Artifact._save_skip_storage = _save_skip_storage
Artifact._cache_path = _cache_path
Artifact.path = path
Artifact.describe = describe
Artifact.view_lineage = view_lineage
