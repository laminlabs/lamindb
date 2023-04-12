import traceback
from functools import partial
from typing import List, Optional, Tuple, Union, overload  # noqa

import lndb
import sqlmodel as sqm
from lamin_logger import logger
from lndb import settings as setup_settings
from lndb_storage import delete_storage, store_object, write_adata_zarr
from lndb_storage._file import print_hook
from lnschema_core import File
from lnschema_core.dev._storage import storage_key_from_file
from pydantic.fields import ModelPrivateAttr

from .._docs import doc_args
from ._core import file_to_sqm, get_session_from_kwargs
from ._select import select

add_docs = """
Insert or update data records.

Inserts a new :term:`record` if the corresponding row doesn't exist.
Updates the corresponding row with the record if it exists.

To update a row, query it with `.select` and modify it before
passing it to `add`.

Guide: :doc:`/guide/add-delete`.

Args:
    record: One or multiple records as instances of `SQLModel`.

Returns:
    The record as returned from the database with a `created_at` timestamp.

Examples:

    Add a record (errors if already exists):

    >>> ln.add(ln.Transform(name="My pipeline"))
    Transform(id="0Cb86EZj", name="My pipeline", ...)

    Update an existing record:

    >>> transform = ln.select(ln.Transform, id="0Cb86EZj").one()
    >>> transform.name = "New name"
    >>> ln.add(transform)
    Transform(id="0Cb86EZj", name="New name", ...)

    Add a record with passed fields if not yet exists:

    >>> # add a record if the metadata combination is not already exist in the DB
    >>> # if exists, returns the existing record from the DB
    >>> ln.add(ln.Transform, name="My transform", v="1")
    Transform(id="0Cb86EZj", name="My pipeline", ...)
    >>> # is equivalent to the following:
    >>> transform = ln.select(ln.Transform, name="My transform", v="1").one_or_none()
    >>> if transform is None:
    >>>     ln.add(transform)

"""


@overload
def add(record: sqm.SQLModel) -> sqm.SQLModel:
    ...


# Currently seeing the following error without type ignore:
# Overloaded function signature 2 will never be matched: signature 1's parameter
# type(s) are the same or broader
@overload
def add(records: List[sqm.SQLModel]) -> List[sqm.SQLModel]:  # type: ignore
    ...


@overload
def add(  # type: ignore
    entity: sqm.SQLModel, **fields
) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
    ...


@doc_args(add_docs)
def add(  # type: ignore
    record: Union[sqm.SQLModel, List[sqm.SQLModel]], **fields
) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
    """{}"""  # noqa
    session = get_session_from_kwargs(fields)
    if isinstance(record, list):
        records = record
    elif isinstance(record, sqm.SQLModel):
        records = [record]
    else:
        model = file_to_sqm(record)
        results = select(model, **fields).one_or_none()
        if results is None:
            records = [model(**fields)]
        else:
            logger.info(
                f"An existing {results.__class__.__name__} record is found in the DB:"
            )
            return results

    if session is None:  # assume global session
        session = setup_settings.instance.session()
        setup_settings.instance._cloud_sqlite_locker.lock()
        close = True
    else:
        close = False

    # commit metadata to database
    db_error = None
    for record in records:
        # the following ensures that queried objects (within __init__)
        # behave like queried objects, only example right now: Run
        if (
            _private_not_empty(record, "_ln_identity_key")
            and record._ln_identity_key is not None  # noqa
        ):
            record._sa_instance_state.key = record._ln_identity_key
        session.add(record)
    try:
        session.commit()
    except Exception as e:
        db_error = e

    # upload data objects to storage
    added_records = []
    if db_error is None:
        added_records, upload_error = upload_committed_records(records, session)

    if close:
        session.close()
        setup_settings.instance._update_cloud_sqlite_file()
        setup_settings.instance._cloud_sqlite_locker.unlock()

    error = db_error or upload_error
    if error is not None:
        error_message = prepare_error_message(records, added_records, error)
        raise RuntimeError(error_message)
    elif len(added_records) > 1:
        return added_records
    else:
        return added_records[0]


def upload_committed_records(records, session):
    """Upload records in a list of database-commited records to storage.

    If any upload fails, subsequent records are cleaned up from the DB.
    """
    # make sure ALL records are up-to-date to enable accurate comparison
    # during metadata cleanup
    error = None
    added_records = []
    for record in records:
        session.refresh(record)

    # upload data objects
    for record in records:
        if isinstance(record, File) and _private_not_empty(record, "_local_filepath"):
            try:
                upload_data_object(record)
            except Exception as e:
                error = e
                break
        added_records += [record]

    # clear old files on update
    for record in added_records:
        if isinstance(record, File) and _private_not_empty(record, "_clear_storagekey"):
            try:
                if record._clear_storagekey is not None:
                    delete_storage(record._clear_storagekey)
                    record._clear_storagekey = None
            except Exception as e:
                error = e

    # clean up metadata for objects not uploaded to storage
    if error is not None:
        for record in records:
            if record not in added_records:
                session.delete(record)
        session.commit()

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
    if (
        _private_not_empty(file, "_to_store")
        and file._to_store
        and file.suffix != ".zarr"
    ):
        logger.hint(f"storing object {file.name} with key {file_storage_key}")
        store_object(file._local_filepath, file_storage_key)
    elif (
        file.suffix == ".zarr"
        and _private_not_empty(file, "_memory_rep")
        and file._memory_rep is not None
    ):
        logger.hint(f"storing object {file.name} with key {file_storage_key}")
        storagepath = lndb.settings.storage.key_to_filepath(file_storage_key)
        print_progress = partial(print_hook, filepath=file.name)
        write_adata_zarr(file._memory_rep, storagepath, callback=print_progress)


# sometimes obj.attr results in
def _private_not_empty(obj, attr):
    if hasattr(obj, attr):
        return not isinstance(getattr(obj, attr), ModelPrivateAttr)
    else:
        return False
