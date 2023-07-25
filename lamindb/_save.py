import os
import shutil
import traceback
from datetime import datetime
from functools import partial
from typing import Iterable, List, Optional, Tuple, Union, overload  # noqa

import lamindb_setup
from django.db import transaction
from lamin_utils import logger
from lnschema_core.models import ORM, File

from lamindb.dev.storage import store_object
from lamindb.dev.storage.file import (
    auto_storage_key_from_file,
    delete_storage_using_key,
    print_hook,
)

try:
    from lamindb.dev.storage._zarr import write_adata_zarr
except ImportError:

    def write_adata_zarr(filepath):  # type: ignore
        raise ImportError("Please install zarr: pip install zarr")


def save(records: Iterable[ORM], **kwargs) -> None:  # type: ignore
    """Bulk save to database & storage.

    Note:

        This is a much **faster** way to save many records in the database.

    Warning:

        It neither automatically creates related records nor updates existing records!
        Use ``ORM.save()`` for these use cases.

    Args:
        records: One or multiple ``ORM`` objects.

    Examples:

        Save a collection of records in one transaction, which is much faster
        than writing a loop over calls ``projects.save()``:

        >>> labels = [ln.Label(f"Label {i}") for i in range(10)]
        >>> ln.save(projects)

        For a single record, use ``.save()``:

        >>> transform = ln.Transform(name="My pipeline")
        >>> transform.save()

        Update a single existing record:

        >>> transform = ln.select(ln.Transform, id="0Cb86EZj").one()
        >>> transform.name = "New name"
        >>> transform.save()

    """
    if isinstance(records, Iterable):
        records = set(records)
    elif isinstance(records, ORM):
        records = {records}

    # we're distinguishing between files and non-files
    # because for files, we want to bulk-upload
    # rather than upload one-by-one
    files = {r for r in records if isinstance(r, File)}
    non_files = records.difference(files)
    if non_files:
        # first save all records without recursing parents
        bulk_create(non_files)
        non_files_with_parents = {r for r in non_files if hasattr(r, "_parents")}

        if len(non_files_with_parents) > 0 and kwargs.get("parents") is not False:
            # this can only happen within lnschema_bionty right now!!
            # we might extend to core lamindb later
            import lnschema_bionty as lb

            if kwargs.get("parents") or (
                kwargs.get("parents") is None and lb.settings.auto_save_parents
            ):
                # save the record with parents one by one
                logger.warning(
                    "Now recursing through parents: "
                    "this only happens once, but is much slower than bulk saving"
                )
                logger.hint(
                    "You can switch this off via: lb.settings.auto_save_parents = False"
                )
                for record in non_files_with_parents:
                    record._save_ontology_parents()

    if files:
        with transaction.atomic():
            for record in files:
                record._save_skip_storage()
        store_files(files)

    # this function returns None as potentially 10k records might be saved
    # refreshing all of them from the DB would mean a severe performance penalty
    # 2nd reason: consistency with Django Model.save(), which also returns None
    return None


def bulk_create(records: Iterable[ORM]):
    orm = next(iter(records)).__class__
    orm.objects.bulk_create(records, ignore_conflicts=True)


# This is also used within File.save()
def check_and_attempt_upload(file: File) -> Optional[Exception]:
    # if File object is either newly instantiated or replace() was called on
    # a local env it will have a _local_filepath and needs to be uploaded
    if hasattr(file, "_local_filepath"):
        try:
            upload_data_object(file)
        except Exception as exception:
            logger.warning(f"Could not upload file: {file}")
            return exception
        # copies (if ob-disk) or moves the temporary file (if in-memory) to the cache
        copy_or_move_to_cache(file)
        # after successful upload, we should remove the attribute so that another call
        # call to save won't upload again, the user should call replace() then
        del file._local_filepath
    # returning None means proceed (either success or no action needed)
    return None


def copy_or_move_to_cache(file: File):
    local_path = file._local_filepath

    # in-memory zarr or not cloud
    if local_path is None or not lamindb_setup.settings.storage.is_cloud:
        return None

    local_path = local_path.resolve()

    # on-disk zarr
    if not local_path.is_file():
        return None

    # maybe create something like storage.key_to_local(key) later to simplfy
    storage_key = auto_storage_key_from_file(file)
    storage_path = lamindb_setup.settings.storage.key_to_filepath(storage_key)
    cache_path = lamindb_setup.settings.storage.cloud_to_local_no_update(storage_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    cache_dir = lamindb_setup.settings.storage.cache_dir

    if cache_dir in local_path.parents:
        local_path.replace(cache_path)
    else:
        shutil.copy(local_path, cache_path)
    # make sure that the cached version is older than the cloud one
    mts = datetime.now().timestamp() + 1.0
    os.utime(cache_path, times=(mts, mts))


# This is also used within File.save()
def check_and_attempt_clearing(file: File) -> Optional[Exception]:
    # this is a clean-up operation after replace() was called
    # this will only evaluate to True if replace() was called
    if hasattr(file, "_clear_storagekey"):
        try:
            if file._clear_storagekey is not None:
                delete_storage_using_key(file, file._clear_storagekey)
                logger.success(
                    f"Deleted stale object at storage key {file._clear_storagekey}"
                )
                file._clear_storagekey = None
        except Exception as exception:
            return exception
    # returning None means proceed (either success or no action needed)
    return None


def store_files(files: Iterable[File]) -> None:
    """Upload files in a list of database-committed files to storage.

    If any upload fails, subsequent files are cleaned up from the DB.
    """
    exception: Optional[Exception] = None
    # because uploads might fail, we need to maintain a new list
    # of the succeeded uploads
    stored_files = []

    # upload new local files
    for file in files:
        exception = check_and_attempt_upload(file)
        if exception is not None:
            break
        stored_files += [file]
        exception = check_and_attempt_clearing(file)
        if exception is not None:
            logger.warning(f"clean up of {file._clear_storagekey} failed")
            break

    if exception is not None:
        # clean up metadata for files not uploaded to storage
        with transaction.atomic():
            for file in files:
                if file not in stored_files:
                    file._delete_skip_storage()
        error_message = prepare_error_message(files, stored_files, exception)
        raise RuntimeError(error_message)
    return None


def prepare_error_message(records, stored_files, exception) -> str:
    if len(records) == 1 or len(stored_files) == 0:
        error_message = (
            "No entries were uploaded or committed"
            " to the database. See error message:\n\n"
        )
    else:
        error_message = (
            "The following entries have been"
            " successfully uploaded and committed to the database:\n"
        )
        for record in stored_files:
            error_message += (
                f"- {', '.join(record.__repr__().split(', ')[:3]) + ', ...)'}\n"
            )
        error_message += "\nSee error message:\n\n"
    error_message += f"{str(exception)}\n\n{traceback.format_exc()}"
    return error_message


def upload_data_object(file) -> None:
    """Store and add file and its linked entries."""
    # do NOT hand-craft the storage key!
    file_storage_key = auto_storage_key_from_file(file)
    if hasattr(file, "_to_store") and file._to_store and file.suffix != ".zarr":
        logger.hint(f"storing file {file.id} with key {file_storage_key}")
        store_object(file._local_filepath, file_storage_key)
    elif (
        file.suffix in {".zarr", ".zrad"}
        and hasattr(file, "_memory_rep")
        and file._memory_rep is not None
    ):
        logger.hint(f"storing file {file.id} with key {file_storage_key}")
        storagepath = lamindb_setup.settings.storage.key_to_filepath(file_storage_key)
        print_progress = partial(print_hook, filepath=file_storage_key)
        write_adata_zarr(file._memory_rep, storagepath, callback=print_progress)
