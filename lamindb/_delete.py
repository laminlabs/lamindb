import sqlmodel as sqm
from lndb_setup import settings
from lnschema_core import DObject, RunIn, Usage

from ._logger import colors, logger
from .dev._core import storage_key_from_dobject
from .dev.file import delete_storage


def delete(record: sqm.SQLModel):
    """Delete data records & data objects.

    Guide: :doc:`/db/guide/add-delete`.

    Example:

    >>> # Delete metadata records
    >>> experiment = ln.select(Experiment, id=experiment_id)
    >>> db.delete(experiment)
    >>> # Delete data objects
    >>> dobject = ln.select(DObject, id=dobject_id)
    >>> db.delete(dobject)

    Args:
        record: One or multiple records as instances of `SQLModel`.
    """
    session = settings.instance.session()
    if isinstance(record, DObject):
        # delete usage events related to the dobject that's to be deleted
        events = session.exec(sqm.select(Usage).where(Usage.dobject_id == record.id))
        for event in events:
            session.delete(event)
        session.commit()
        # delete run_ins related to the dobject that's to be deleted
        run_ins = session.exec(sqm.select(RunIn).where(RunIn.dobject_id == record.id))
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
    if settings.instance._session is None:
        session.close()
    if isinstance(record, DObject):
        # TODO: do not track deletes until we come up
        # with a good design that respects integrity
        # track_usage(entry.id, "delete")
        storage_key = storage_key_from_dobject(record)
        delete_storage(storage_key)
        logger.success(
            f"Deleted {colors.yellow(f'object {storage_key}')} from storage."
        )
