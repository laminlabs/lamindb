from functools import partial
from pathlib import Path
from typing import List, Union, overload  # noqa

import sqlmodel as sqm
from cloudpathlib import CloudPath
from lndb import settings as setup_settings
from lnschema_core import DObject
from sqlalchemy.orm.attributes import set_attribute

from .._docs import doc_args
from ..file import store_file, write_adata_zarr
from ..file._file import print_hook
from ._core import dobject_to_sqm, get_session_from_kwargs
from ._select import select

add_docs = """
Insert or update data records.

Inserts a new :term:`record` if the corresponding row doesn't exist.
Updates the corresponding row with the record if it exists.

To update a row, query it with `.get` or `.select` and modify it before
passing it to `add`.

Guide: :doc:`/guide/add-delete`.

Example:

>>> # add a record (by passing a record)
>>> ln.add(wetlab.Experiment(name="My test", biometa_id=test_id))
>>> # update an existing record
>>> experiment = ln.select(wetlab.Experiment, id=experiment_id).one()
>>> experiment.name = "New name"
>>> ln.add(experiment)
>>> # add a record by fields if not yet exists
>>> ln.add(wetlab.Experiment, name="My test", biometa_id=test_id)

Args:
    record: One or multiple records as instances of `SQLModel`.
    use_fsspec: Whether to use fsspec.
"""


@overload
def add(record: sqm.SQLModel, use_fsspec: bool = True) -> sqm.SQLModel:
    ...


# Currently seeing the following error without type ignore:
# Overloaded function signature 2 will never be matched: signature 1's parameter
# type(s) are the same or broader
@overload
def add(  # type: ignore
    records: List[sqm.SQLModel], use_fsspec: bool = True
) -> List[sqm.SQLModel]:
    ...


@overload
def add(  # type: ignore
    entity: sqm.SQLModel, use_fsspec: bool = True, **fields
) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
    ...


@doc_args(add_docs)
def add(  # type: ignore
    record: Union[sqm.SQLModel, List[sqm.SQLModel]], use_fsspec: bool = True, **fields
) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
    """{}"""  # noqa
    session = get_session_from_kwargs(fields)
    if isinstance(record, list):
        records = record
    elif isinstance(record, sqm.SQLModel):
        records = [record]
    else:
        model = dobject_to_sqm(record)
        results = select(model, **fields).all()
        if len(results) == 1:
            return results[0]
        elif len(results) > 1:
            return results
        else:
            records = [model(**fields)]

    if session is None:  # assume global session
        session = setup_settings.instance.session()
        setup_settings.instance._cloud_sqlite_locker.lock()
        close = True
    else:
        close = False

    # commit metadata to database
    db_error = None
    for record in records:
        prepare_filekey_metadata(record)
        session.add(record)
    try:
        session.commit()
    except Exception as e:
        db_error = e

    # upload data objects to storage
    added_records = []
    if db_error is None:
        added_records, upload_error = upload_committed_records(
            records, session, use_fsspec=use_fsspec
        )

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


def upload_committed_records(records, session, use_fsspec):
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
        if isinstance(record, DObject) and hasattr(record, "_local_filepath"):
            try:
                upload_data_object(record, use_fsspec=use_fsspec)
            except Exception as e:
                error = e
                break
        added_records += [record]

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
            "An error occured. No entries were uploaded or committed"
            " to the database. See error message below.\n\n"
        )
    else:
        error_message = (
            "An error occured. The following entries have been"
            " successfully uploaded and committed to the database:\n"
        )
        for record in added_records:
            error_message += (
                f"- {', '.join(record.__repr__().split(', ')[:3]) + ', ...)'}\n"
            )
        error_message += "\nSee error message below.\n\n"
    error_message += f"{str(error)}"
    return error_message


def prepare_filekey_metadata(record) -> None:
    """For cloudpath, write custom filekey to _filekey.

    A filekey excludes the storage root and the file suffix.
    """

    def set_filekey(record: sqm.SQLModel, filepath: Union[Path, CloudPath]):
        set_attribute(
            record,
            "_filekey",
            str(filepath)
            .replace(f"{setup_settings.instance.storage.root}/", "")
            .replace(record.suffix, ""),
        )

    # _local_filepath private attribute is only added
    # when creating DObject from data
    if isinstance(record, DObject) and hasattr(record, "_local_filepath"):
        if record.suffix == ".zarr":
            pass
        else:
            # cloud storage
            if record._cloud_filepath is not None:
                set_filekey(record, record._cloud_filepath)
            # local storage that is configured
            else:
                if (
                    setup_settings.instance.storage.root
                    in record._local_filepath.parents
                ):
                    set_filekey(record, record._local_filepath)


def upload_data_object(dobject, use_fsspec: bool = True) -> None:
    """Store and add dobject and its linked entries."""
    dobject_storage_key = f"{dobject.id}{dobject.suffix}"

    if dobject.suffix != ".zarr":
        # no file upload for cloud storage or configured local storage
        # only copy data when file is not in the configured local storage
        if (dobject._cloud_filepath is None) and (
            setup_settings.instance.storage.root not in dobject._local_filepath.parents
        ):
            store_file(
                dobject._local_filepath, dobject_storage_key, use_fsspec=use_fsspec
            )
    else:
        storagepath = setup_settings.instance.storage.key_to_filepath(
            dobject_storage_key
        )
        print_progress = partial(print_hook, filepath=dobject.name)
        write_adata_zarr(dobject._memory_rep, storagepath, callback=print_progress)
