from lamin_logger import logger
from lamindb_setup import settings
from lamindb_setup.dev import StorageSettings
from lnschema_core.models import File, Storage


# add type annotations back asap when re-organizing the module
def auto_storage_key_from_file(file: File):
    if file.key is None:
        return f"lndb/{file.id}{file.suffix}"
    else:
        return file.key


def attempt_accessing_path(file: File, storage_key: str):
    if file.storage_id == settings.storage.id:
        path = settings.storage.key_to_filepath(storage_key)
    else:
        logger.warning(
            "file.path() is slow for files outside the currently configured storage"
            " location\nconsider joining for the set of files you're interested in:"
            " ln.select(ln.File, ln.Storage)the path is storage.root / file.key if"
            " file.key is not None\notherwise storage.root / (file.id + file.suffix)"
        )
        storage = Storage.select(id=file.storage_id).one()
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
        return attempt_accessing_path(file, legacy_storage_key)
    return path


# add type annotations back asap when re-organizing the module
def filepath_from_file(file: File):
    storage_key = auto_storage_key_from_file(file)
    path = attempt_accessing_path(file, storage_key)
    return path
