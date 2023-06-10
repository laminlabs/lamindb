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

    added_records, upload_error = store_files(records)

    if upload_error is not None:
        error_message = prepare_error_message(records, added_records, upload_error)
        raise RuntimeError(error_message)

    # this function returns None as potentially 10k records might be saved
    # refreshing all of them from the DB would mean a severe performance penatly
    # 2nd reason: consistency with Django Model.save(), which also returns None


def store_files(records):
    """Upload records in a list of database-committed records to storage.

    If any upload fails, subsequent records are cleaned up from the DB.
    """
    error = None
    added_records = []

    # upload data objects
    for record in records:
        if isinstance(record, File) and hasattr(record, "_local_filepath"):
            try:
                upload_data_object(record)
            except Exception as e:
                logger.warning(f"Could not upload file: {record}")
                error = e
                break
        added_records += [record]

    # clear old files on update
    for record in added_records:
        if isinstance(record, File) and hasattr(record, "_clear_storagekey"):
            try:
                if record._clear_storagekey is not None:
                    delete_storage(record._clear_storagekey)
                    record._clear_storagekey = None
            except Exception as e:
                error = e

    # clean up metadata for objects not uploaded to storage
    if error is not None:
        with transaction.atomic():
            for record in records:
                if record not in added_records:
                    record.delete()

    return added_records, error


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
