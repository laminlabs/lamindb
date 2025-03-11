from __future__ import annotations

import shutil
from typing import TYPE_CHECKING

import fsspec
from lamindb_setup.core import StorageSettings
from lamindb_setup.core.upath import (
    LocalPathClasses,
    UPath,
)

from lamindb.core._settings import settings

if TYPE_CHECKING:
    from pathlib import Path

    from lamindb_setup.core.types import UPathStr

    from lamindb.models.artifact import Artifact


AUTO_KEY_PREFIX = ".lamindb/"


# add type annotations back asap when re-organizing the module
def auto_storage_key_from_artifact(artifact: Artifact):
    if artifact.key is None or artifact._key_is_virtual:
        is_dir = artifact.n_files is not None
        return auto_storage_key_from_artifact_uid(artifact.uid, artifact.suffix, is_dir)
    else:
        return artifact.key


def auto_storage_key_from_artifact_uid(uid: str, suffix: str, is_dir: bool) -> str:
    assert isinstance(suffix, str)  # noqa: S101 Suffix cannot be None.
    if is_dir:
        uid_storage = uid[:16]  # 16 chars, leave 4 chars for versioning
    else:
        uid_storage = uid
    storage_key = f"{AUTO_KEY_PREFIX}{uid_storage}{suffix}"
    return storage_key


def check_path_is_child_of_root(path: UPathStr, root: UPathStr) -> bool:
    if fsspec.utils.get_protocol(str(path)) != fsspec.utils.get_protocol(str(root)):
        return False
    path_upath = UPath(path)
    root_upath = UPath(root)
    if path_upath.protocol == "s3":
        endpoint_path = path_upath.storage_options.get("endpoint_url", "")
        endpoint_root = root_upath.storage_options.get("endpoint_url", "")
        if endpoint_path != endpoint_root:
            return False
    # we don't resolve http links because they can resolve into a different domain
    # for example into a temporary url
    if path_upath.protocol not in {"http", "https"}:
        path_upath = path_upath.resolve()
        root_upath = root_upath.resolve()
    # str is needed to eliminate UPath storage_options
    # which affect equality checks
    return UPath(str(root_upath)) in UPath(str(path_upath)).parents


# returns filepath and root of the storage
def attempt_accessing_path(
    artifact: Artifact,
    storage_key: str,
    using_key: str | None = None,
    access_token: str | None = None,
) -> tuple[UPath, StorageSettings]:
    # check whether the file is in the default db and whether storage
    # matches default storage
    from lamindb.models import Storage

    if (
        artifact._state.db in ("default", None)
        and artifact.storage_id == settings._storage_settings.id
    ):
        if access_token is None:
            storage_settings = settings._storage_settings
        else:
            storage_settings = StorageSettings(
                settings.storage.root, access_token=access_token
            )
    else:
        if artifact._state.db not in ("default", None) and using_key is None:
            storage = Storage.using(artifact._state.db).get(id=artifact.storage_id)
        else:
            storage = Storage.objects.using(using_key).get(id=artifact.storage_id)
        # find a better way than passing None to instance_settings in the future!
        storage_settings = StorageSettings(storage.root, access_token=access_token)
    path = storage_settings.key_to_filepath(storage_key)
    return path, storage_settings


def filepath_from_artifact(
    artifact: Artifact, using_key: str | None = None
) -> tuple[UPath, StorageSettings | None]:
    if hasattr(artifact, "_local_filepath") and artifact._local_filepath is not None:
        return artifact._local_filepath.resolve(), None
    storage_key = auto_storage_key_from_artifact(artifact)
    path, storage_settings = attempt_accessing_path(
        artifact, storage_key, using_key=using_key
    )
    return path, storage_settings


# virtual key is taken into consideration
# only if the version is latest
def _cache_key_from_artifact_storage(
    artifact: Artifact, storage_settings: StorageSettings | None
):
    cache_key = None
    if (
        artifact._key_is_virtual
        and artifact.key is not None
        and storage_settings is not None
        and artifact.is_latest
    ):
        cache_key = (storage_settings.root / artifact.key).path
    return cache_key


# return filepath and cache_key if needed
def filepath_cache_key_from_artifact(
    artifact: Artifact, using_key: str | None = None
) -> tuple[UPath, str | None]:
    filepath, storage_settings = filepath_from_artifact(artifact, using_key)
    if isinstance(filepath, LocalPathClasses):
        return filepath, None
    cache_key = _cache_key_from_artifact_storage(artifact, storage_settings)
    return filepath, cache_key


def store_file_or_folder(
    local_path: UPathStr, storage_path: UPath, print_progress: bool = True, **kwargs
) -> None:
    """Store file or folder (localpath) at storagepath."""
    local_path = UPath(local_path)
    if not isinstance(storage_path, LocalPathClasses):
        # this uploads files and directories
        if local_path.is_dir():
            create_folder = False
            try:
                # if storage_path already exists we need to delete it
                # if local_path is a directory
                # to replace storage_path correctly
                if storage_path.stat().as_info()["type"] == "directory":
                    storage_path.rmdir()
                else:
                    storage_path.unlink()
            except (FileNotFoundError, PermissionError):
                pass
        else:
            create_folder = None
        storage_path.upload_from(
            local_path,
            create_folder=create_folder,
            print_progress=print_progress,
            **kwargs,
        )
    else:  # storage path is local
        if local_path.resolve().as_posix() == storage_path.resolve().as_posix():
            return None
        storage_path.parent.mkdir(parents=True, exist_ok=True)
        if local_path.is_file():
            shutil.copyfile(local_path, storage_path)
        else:
            if storage_path.exists():
                shutil.rmtree(storage_path)
            shutil.copytree(local_path, storage_path)


def delete_storage_using_key(
    artifact: Artifact,
    storage_key: str,
    raise_file_not_found_error: bool = True,
    using_key: str | None = None,
) -> None | str:
    filepath, _ = attempt_accessing_path(artifact, storage_key, using_key=using_key)
    return delete_storage(
        filepath, raise_file_not_found_error=raise_file_not_found_error
    )


def delete_storage(
    storagepath: Path, raise_file_not_found_error: bool = True
) -> None | str:
    """Delete arbitrary artifact."""
    if storagepath.is_file():
        storagepath.unlink()
    elif storagepath.is_dir():
        if isinstance(storagepath, LocalPathClasses) or not isinstance(
            storagepath, UPath
        ):
            shutil.rmtree(storagepath)
        else:
            storagepath.rmdir()
    elif raise_file_not_found_error:
        raise FileNotFoundError(f"{storagepath} is not an existing path!")
    else:
        return "did-not-delete"
    return None
