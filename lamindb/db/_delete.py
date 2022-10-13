import lnschema_core as schema_core
import sqlmodel as sqm
from lndb_setup import settings

from .._logger import colors, logger
from ..dev._core import storage_key_from_dobject
from ..dev.file import delete_file
from ..schema._table import Table


def _create_delete_func(model):
    name = model.__name__

    def delete_func(key):
        with sqm.Session(settings.instance.db_engine()) as session:
            entry = session.get(model, key)
            if entry is None:
                logger.warning(f"Entry {key} does not exist.")
                return None
            if name == "dobject":
                # delete usage events related to the dobject that's to be deleted
                events = session.exec(
                    sqm.select(schema_core.usage).where(
                        schema_core.usage.dobject_id == key
                    )
                )
                for event in events:
                    session.delete(event)
                session.commit()
                # delete dtransform_ins related to the dobject that's to be deleted
                dtransform_ins = session.exec(
                    sqm.select(schema_core.dtransform_in).where(
                        schema_core.dtransform_in.dobject_id == key
                    )
                )
                for dtransform_in in dtransform_ins:
                    session.delete(dtransform_in)
                session.commit()
            session.delete(entry)
            session.commit()
            settings.instance._update_cloud_sqlite_file()
            logger.success(
                f"Deleted {colors.yellow(f'entry {key}')} in"
                f" {colors.blue(f'table {name}')}."
            )
        if name == "dobject":
            # TODO: do not track deletes until we come up
            # with a good design that respects integrity
            # track_usage(entry.id, "delete")
            storage_key = storage_key_from_dobject(entry)
            delete_file(storage_key)
            logger.success(
                f"Deleted {colors.yellow(f'object {storage_key}')} from storage."
            )

    delete_func.__name__ = name
    return delete_func


class delete:
    """Delete data.

    Deletes the row identified by the key in the table that represents the entity.
    For `dobject`, this also deletes the data object in storage.

    Example:

    >>> delete.dobject("MHnsI0TWWNu3VBsedg6gS")
    """

    pass


for model in Table.list_models():
    func = _create_delete_func(model=model)
    setattr(delete, model.__name__, staticmethod(func))
