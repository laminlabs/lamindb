from typing import Union

from lamin_logger import logger
from lamindb_setup import settings
from lamindb_setup.dev import StorageSettings
from lnschema_core.models import File, Folder, Storage


# add type annotations back asap when re-organizing the module
def storage_key_from_file(file: File):
    if file.key is None:
        return f"lndb/{file.id}{file.suffix}"
    else:
        return file.key


def attempt_accessing_path(file_or_folder: Union[File, Folder], storage_key: str):
    if file_or_folder.storage_id == settings.storage.id:
        path = settings.storage.key_to_filepath(storage_key)
    else:
        logger.warning(
            "file.path() is slow for files outside the currently configured storage"
            " location\nconsider joining for the set of files you're interested in:"
            " ln.select(ln.File, ln.Storage)the path is storage.root / file.key if"
            " file.key is not None\notherwise storage.root / (file.id + file.suffix)"
        )
        storage = Storage.select(id=file_or_folder.storage_id).one()
        # find a better way than passing None to instance_settings in the future!
        storage_settings = StorageSettings(storage.root, instance_settings=None)
        path = storage_settings.key_to_filepath(storage_key)
    # the following is for backward compat
    if storage_key.startswith("lndb/") and not path.exists():
        logger.warning(
            "You have auto-keyed files in your storage root, please move them into an"
            " 'lndb/' subfolder"
        )
        legacy_storage_key = storage_key.lstrip("/lndb")
        attempt_accessing_path(file_or_folder, legacy_storage_key)
    return path


# add type annotations back asap when re-organizing the module
def filepath_from_file_or_folder(file_or_folder: Union[File, Folder]):
    if isinstance(file_or_folder, File):
        storage_key = storage_key_from_file(file_or_folder)
    else:
        storage_key = file_or_folder.key
        if storage_key is None:
            raise ValueError("Only real folders have a path!")
    path = attempt_accessing_path(file_or_folder, storage_key)
    return path
