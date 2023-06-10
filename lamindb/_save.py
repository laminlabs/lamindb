import traceback
from functools import partial
from typing import Iterable, List, Optional, Tuple, Union, overload  # noqa

import lamindb_setup
from django.db import transaction
from lamin_logger import logger
from lnschema_core.models import BaseORM, File, storage_key_from_file

from lamindb.dev.storage import delete_storage, store_object, write_adata_zarr
from lamindb.dev.storage._file import print_hook


@overload
def save(record: BaseORM) -> BaseORM:
    ...


# Currently seeing the following error without type ignore:
# Overloaded function signature 2 will never be matched: signature 1's parameter
# type(s) are the same or broader
@overload
def save(records: Iterable[BaseORM]) -> Iterable[BaseORM]:  # type: ignore
    ...


def save(record: Union[BaseORM, Iterable[BaseORM]], **fields) -> None:  # type: ignore
    """Save to database & storage.

    Inserts a new :term:`record` if the corresponding row doesn't exist.
    Updates the corresponding row with the record if it exists.

    To update a row, query it with `.select` and modify it before
    passing it to `save`.

    Args:
        record: One or multiple `BaseORM` objects.

    Returns:
        The record as returned from the database with a `created_at` timestamp.

    Examples:

        Save a record (errors if already exists):

        >>> transform = ln.Transform(name="My pipeline")
        >>> ln.save(transform)  # equivalent to transform.save()

        Save a collection of records in one transaction:

        >>> projects = [ln.Project(f"Project {i}") for i in range(10)]
        >>> ln.save(projects)

        Update an existing record:

        >>> transform = ln.select(ln.Transform, id="0Cb86EZj").one()
        >>> transform.name = "New name"
        >>> ln.save(transform)  # equivalent to transform.save()

    """
    if isinstance(record, Iterable):
        records = set(record)
    elif isinstance(record, BaseORM):
        records = {record}

    # we're distinguishing between files and non-files
    # because for files, we want to bulk-upload
    # rather than upload one-by-one
    files = {r for r in records if isinstance(r, File)}
    non_files = records.difference(files)
    if non_files:
        with transaction.atomic():
            for record in non_files:
                record.save()
    if files:
        with transaction.atomic():
            for record in files:
                record._save_kip_store()
        store_files(files)

    # this function returns None as potentially 10k records might be saved
    # refreshing all of them from the DB would mean a severe performance penalty
    # 2nd reason: consistency with Django Model.save(), which also returns None
    return None


def check_and_attempt_upload(file: File) -> Optional[Exception]:
    # if File object is either newly instantiated or replace() was called on
    # a local env it will have a _local_filepath and needs to be uploaded
    if hasattr(file, "_local_filepath"):
        try:
            upload_data_object(file)
        except Exception as exception:
            logger.warning(f"Could not upload file: {file}")
            return exception
    # returning None means proceed (either success or no action needed)
    return None


def check_and_attempt_clearing(file: File) -> Optional[Exception]:
    # this is a clean-up operation after replace() was called
    # this will only evaluate to True if replace() was called
    if hasattr(file, "_clear_storagekey"):
        try:
            if file._clear_storagekey is not None:
                delete_storage(file._clear_storagekey)
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
                    file.delete()
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
    file_storage_key = storage_key_from_file(file)
    if hasattr(file, "_to_store") and file._to_store and file.suffix != ".zarr":
        logger.hint(f"storing object {file.name} with key {file_storage_key}")
        store_object(file._local_filepath, file_storage_key)
    elif (
        file.suffix == ".zarr"
        and hasattr(file, "_memory_rep")
        and file._memory_rep is not None
    ):
        logger.hint(f"storing object {file.name} with key {file_storage_key}")
        storagepath = lamindb_setup.settings.storage.key_to_filepath(file_storage_key)
        print_progress = partial(print_hook, filepath=file.name)
        write_adata_zarr(file._memory_rep, storagepath, callback=print_progress)
