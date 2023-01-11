from typing import List, Union, overload  # noqa

import sqlmodel as sqm
from lndb_setup import settings
from lnschema_core import DObject, Usage
from lnschema_core.link import RunIn

from ._logger import colors, logger
from .dev._core import storage_key_from_dobject
from .dev.db._select import select
from .dev.file import delete_storage


@overload
def delete(record: sqm.SQLModel) -> None:
    ...


@overload
def delete(records: List[sqm.SQLModel]) -> None:  # type: ignore
    ...


@overload
def delete(entity: sqm.SQLModel, **fields) -> None:  # type: ignore
    ...


def delete(  # type: ignore
    record: Union[sqm.SQLModel, List[sqm.SQLModel]], **fields
) -> None:
    """Delete data records & data objects.

    Guide: :doc:`/guide/add-delete`.

    Example:

    >>> # Delete by record
    >>> experiment = ln.select(Experiment, id=experiment_id)
    >>> ln.delete(experiment)
    >>> # Delete data objects
    >>> dobject = ln.select(DObject, id=dobject_id)
    >>> ln.delete(dobject)
    >>> # Delete by fields
    >>> ln.delete(DObject, id=dobject_id)

    Args:
        record: One or multiple records as instances of `SQLModel`.
    """
    if isinstance(record, list):
        records = record
    elif isinstance(record, sqm.SQLModel):
        records = [record]
    else:
        model = record
        results = select(model, **fields).all()
        if len(results) == 0:
            return None
        records = results
    session = settings.instance.session()
    for record in records:
        if isinstance(record, DObject):
            # delete usage events related to the dobject that's to be deleted
            events = session.exec(
                sqm.select(Usage).where(Usage.dobject_id == record.id)
            )
            for event in events:
                session.delete(event)
            session.commit()
            # delete run_ins related to the dobject that's to be deleted
            run_ins = session.exec(
                sqm.select(RunIn).where(RunIn.dobject_id == record.id)
            )
            for run_in in run_ins:
                session.delete(run_in)
            session.commit()
        session.delete(record)
        session.commit()
        settings.instance._update_cloud_sqlite_file()
        logger.success(
            f"Deleted {colors.yellow(f'row {record}')} in"
            f" {colors.blue(f'table {type(record).__name__}')}."
        )
        if isinstance(record, DObject):
            # TODO: do not track deletes until we come up
            # with a good design that respects integrity
            # track_usage(entry.id, "delete")
            storage_key = storage_key_from_dobject(record)
            delete_storage(storage_key)
            logger.success(
                f"Deleted {colors.yellow(f'object {storage_key}')} from storage."
            )
    session.close()
