from lamin_logger import logger
from lamindb_setup import settings
from lamindb_setup.dev import StorageSettings
from lnschema_core.models import File, Storage

AUTO_KEY_PREFIX = ".lamindb/"


# add type annotations back asap when re-organizing the module
def auto_storage_key_from_file(file: File):
    if file.key is None:
        return f"{AUTO_KEY_PREFIX}{file.id}{file.suffix}"
    else:
        return file.key


def attempt_accessing_path(file: File, storage_key: str):
    if file.storage_id == settings.storage.id:
        path = settings.storage.key_to_filepath(storage_key)
    else:
        logger.warning(
            "file.path() is slower for files outside the currently configured storage"
            " location"
        )
        storage = Storage.select(id=file.storage_id).one()
        # find a better way than passing None to instance_settings in the future!
        storage_settings = StorageSettings(storage.root, instance_settings=None)
        path = storage_settings.key_to_filepath(storage_key)
    # the following is for backward compat
    if storage_key.startswith(AUTO_KEY_PREFIX) and not path.exists():
        logger.warning(
            "You have auto-keyed files in your storage root, please move them into"
            f" {AUTO_KEY_PREFIX} within your storage location"
        )
        # try legacy_storage_key in root
        for previous_prefix in ["", "lndb/"]:
            legacy_storage_key = storage_key.replace(AUTO_KEY_PREFIX, previous_prefix)
            path = settings.storage.key_to_filepath(legacy_storage_key)
            if path.exists():
                return path
    return path


# add type annotations back asap when re-organizing the module
def filepath_from_file(file: File):
    storage_key = auto_storage_key_from_file(file)
    path = attempt_accessing_path(file, storage_key)
    return path
