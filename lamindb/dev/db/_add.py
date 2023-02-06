from functools import partial
from typing import Dict, List, Tuple, Union, overload  # noqa

import sqlmodel as sqm
from lndb_setup import settings
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
        session = settings.instance.session()
        settings.instance._cloud_sqlite_locker.lock()
        close = True
    else:
        close = False

    # commit metadata to database
    raised_errors = []
    for record in records:
        prepare_filekey_metadata(record)
        session.add(record)
    try:
        session.commit()
    except Exception as e:
        raised_errors.append(e)

    # refresh records and upload data objects to storage
    added_records = []
    if not raised_errors:
        for record in records:
            session.refresh(record)
            if isinstance(record, DObject) and hasattr(record, "_local_filepath"):
                print("Inside LOOP!")
                try:
                    upload_data_object(record, use_fsspec=use_fsspec)
                except Exception as e:
                    # clean up metadata committed to the database
                    print("CAUGHT ERRROR!")
                    raised_errors.append(e)
                    session.delete(record)
                    session.commit()
                    continue
            added_records += [record]

    if close:
        session.close()
        settings.instance._update_cloud_sqlite_file()
        settings.instance._cloud_sqlite_locker.unlock()

    if raised_errors:
        error_message = prepare_error_message(records, added_records)
        raise RuntimeError(error_message)
    elif len(added_records) > 1:
        return added_records
    else:
        return added_records[0]


def prepare_error_message(records, added_records) -> str:
    if len(records) == 1:
        error_message = (
            "An unexpected error occured during upload and no entries were commited to"
            " the database. Please run command again."
        )
    else:
        error_message = (
            "An unexpected error occured during upload.\n\n"
            "The following data objects have been successfully uploaded:\n"
        )
        for record in added_records:
            error_message += (
                f"- {', '.join(record.__repr__().split(', ')[:3]) + ', ...)'}\n"
            )
    return error_message


def prepare_filekey_metadata(record) -> None:
    """For cloudpath, write custom filekey to _filekey."""
    if isinstance(record, DObject) and hasattr(record, "_local_filepath"):
        if record.suffix != ".zarr" and record._cloud_filepath is not None:
            set_attribute(
                record,
                "_filekey",
                str(record._cloud_filepath)
                .replace(f"{settings.instance.storage.root}/", "")
                .split(".")[0],
            )


def upload_data_object(dobject, use_fsspec: bool = True) -> None:
    """Store and add dobject and its linked entries."""
    dobject_storage_key = f"{dobject.id}{dobject.suffix}"

    if dobject.suffix != ".zarr":
        # no file upload for cloud storage
        if dobject._cloud_filepath is None:
            store_file(
                dobject._local_filepath, dobject_storage_key, use_fsspec=use_fsspec
            )
    else:
        storagepath = settings.instance.storage.key_to_filepath(dobject_storage_key)
        print_progress = partial(print_hook, filepath=dobject.name)
        write_adata_zarr(dobject._memory_rep, storagepath, callback=print_progress)
