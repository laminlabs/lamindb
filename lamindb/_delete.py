import traceback
from typing import List, Optional, Union, overload  # noqa

import sqlmodel as sqm
from lndb import settings
from lndb_storage import delete_storage
from lnschema_core import File
from lnschema_core.link import RunIn

from ._logger import colors, logger
from .dev._core import storage_key_from_file
from .dev.db._select import select


@overload
def delete(
    record: sqm.SQLModel,
    delete_data_from_storage: Optional[bool] = None,
) -> None:
    ...


@overload
def delete(
    records: List[sqm.SQLModel],
    delete_data_from_storage: Optional[bool] = None,
) -> None:  # type: ignore
    ...


@overload
def delete(
    entity: sqm.SQLModel,
    delete_data_from_storage: Optional[bool] = None,
    **fields,
) -> None:  # type: ignore
    ...


def delete(  # type: ignore
    record: Union[sqm.SQLModel, List[sqm.SQLModel]],
    delete_data_from_storage: Optional[bool] = None,
    **fields,
) -> None:
    """Delete data records & data objects.

    Guide: :doc:`/guide/add-delete`.

    Example:

    1) Delete by record
    >>> experiment = ln.select(Experiment, id=experiment_id).one()
    >>> ln.delete(experiment)

    2) Delete by fields
    >>> ln.delete(Experiment, id=experiment_id)
    >>> # the result of is equivalent to 1)

    3) Delete data objects (deleting the metadata record and the storage file)
    >>> file = ln.select(File, id=file_id).one()
    >>> # deleting the metadata record occurs automatically
    >>> # you will be asked whether to delete the file from storage
    >>> # or pass boolean values to `delete_data_from_storage`
    >>> ln.delete(file, delete_data_from_storage)

    Args:
        record: One or multiple records as instances of `SQLModel`.
        delete_data_from_storage: Whether to delete data from storage.
    """
    if isinstance(record, list):
        records = record
    elif isinstance(record, sqm.SQLModel):
        records = [record]
    else:
        model = record
        results = select(model, **fields).one_or_none()
        if results is None:
            return results
        else:
            records = [results]

    settings.instance._cloud_sqlite_locker.lock()
    session = settings.instance.session()
    for record in records:
        if isinstance(record, File):
            # delete run_ins related to the file that's to be deleted
            run_ins = session.exec(sqm.select(RunIn).where(RunIn.file_id == record.id))
            for run_in in run_ins:
                session.delete(run_in)
            try:
                session.commit()
            except Exception:
                logger.warning("Deleting run inputs failed!")
                traceback.print_exc()
        session.delete(record)
        try:
            session.commit()
            logger.success(
                f"Deleted {colors.yellow(f'row {record}')} in"
                f" {colors.blue(f'table {type(record).__name__}')}."
            )
        except Exception:
            traceback.print_exc()
        if isinstance(record, File):
            storage_key = storage_key_from_file(record)

            if delete_data_from_storage is None:
                # ask to confirm deleting data from storage
                delete_dialog = (
                    "Confirm Delete: Are you sure you want to delete"
                    f" object {storage_key} from storage? (y/n)"
                )
                decide = input(f"   {delete_dialog}")
            else:
                decide = "y" if delete_data_from_storage else "n"

            if decide not in ("y", "Y", "yes", "Yes", "YES"):
                continue
            try:
                delete_storage(storage_key)
                logger.success(
                    f"Deleted {colors.yellow(f'object {storage_key}')} from storage."
                )
            except Exception:
                traceback.print_exc()
    session.close()
    settings.instance._update_cloud_sqlite_file()
    settings.instance._cloud_sqlite_locker.unlock()
