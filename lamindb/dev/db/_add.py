from functools import partial
from typing import Dict, List, Tuple, Union, overload  # noqa

import sqlmodel as sqm
from lndb_setup import settings
from lnschema_core import DObject

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
    for record in records:
        if isinstance(record, DObject) and hasattr(record, "_local_filepath"):
            upload_data_object(record, use_fsspec=use_fsspec)
    if session is None:  # assume global session
        session = settings.instance.session()
        close = True
    else:
        close = False
    for record in records:
        session.add(record)
    session.commit()
    for record in records:
        session.refresh(record)
    if close:
        session.close()
    settings.instance._update_cloud_sqlite_file()
    if len(records) > 1:
        return records
    else:
        return records[0]


def upload_data_object(dobject, use_fsspec: bool = True) -> None:
    """Store and add dobject and its linked entries."""
    dobject_storage_key = f"{dobject.id}{dobject.suffix}"

    if dobject.suffix != ".zarr":
        store_file(dobject._local_filepath, dobject_storage_key, use_fsspec=use_fsspec)
    else:
        storagepath = settings.instance.storage.key_to_filepath(dobject_storage_key)
        print_progress = partial(print_hook, filepath=dobject._local_filepath)
        write_adata_zarr(dobject._memory_rep, storagepath, callback=print_progress)
