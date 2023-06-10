import traceback
from functools import partial
from typing import List, Optional, Tuple, Union, overload  # noqa

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
def save(records: List[BaseORM]) -> List[BaseORM]:  # type: ignore
    ...


def save(record: Union[BaseORM, List[BaseORM]], **fields) -> None:  # type: ignore
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

        >>> ln.save(ln.Transform(name="My pipeline"))

        Update an existing record:

        >>> transform = ln.select(ln.Transform, id="0Cb86EZj").one()
        >>> transform.name = "New name"
        >>> ln.save(transform)

    """
    if isinstance(record, list):
        records = record
    elif isinstance(record, BaseORM):
        records = [record]

    # no matter which record type, let's commit all of them to the database
    with transaction.atomic():
        for record in records:
            record.save()

    files = [r for r in records if isinstance(r, File)]
    if files:
        store_files(files)

    # this function returns None as potentially 10k records might be saved
    # refreshing all of them from the DB would mean a severe performance penalty
    # 2nd reason: consistency with Django Model.save(), which also returns None


def store_files(files: List[File]):
    """Upload files in a list of database-committed files to storage.

    If any upload fails, subsequent files are cleaned up from the DB.
    """
    error = None
    # because uploads might fail, we need to maintain a new list
    # of the succeeded uploads
    stored_files = []

    # upload new local files
    for file in files:
        # if File object is newly instantiated on a local env
        # it will have a _local_filepath and needs to be uploaded
        if hasattr(file, "_local_filepath"):
            try:
                upload_data_object(file)
            except Exception as e:
                logger.warning(f"Could not upload file: {file}")
                error = e
                break
        # all other files are already stored anyway
        stored_files += [file]

    # clear old files on update
    for record in stored_files:
        if hasattr(record, "_clear_storagekey"):
            try:
                if record._clear_storagekey is not None:
                    delete_storage(record._clear_storagekey)
                    record._clear_storagekey = None
            except Exception as e:
                error = e

    # clean up metadata for objects not uploaded to storage
    if error is not None:
        with transaction.atomic():
            for record in files:
                if record not in stored_files:
                    record.delete()

    error_message = prepare_error_message(files, stored_files, error)
    raise RuntimeError(error_message)


def prepare_error_message(records, added_records, error) -> str:
    if len(records) == 1 or len(added_records) == 0:
        error_message = (
            "No entries were uploaded or committed"
            " to the database. See error message:\n\n"
        )
    else:
        error_message = (
            "The following entries have been"
            " successfully uploaded and committed to the database:\n"
        )
        for record in added_records:
            error_message += (
                f"- {', '.join(record.__repr__().split(', ')[:3]) + ', ...)'}\n"
            )
        error_message += "\nSee error message:\n\n"
    error_message += f"{str(error)}\n\n{traceback.format_exc()}"
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
