from lndb_setup import settings
from sqlmodel import Session

from .._logger import colors, logger
from ..dev import storage_key_from_dobject, track_usage
from ..dev.file import delete_file
from ..schema._table import Table


def _create_delete_func(model):
    name = model.__name__

    def delete_func(key):
        with Session(settings.instance.db_engine()) as session:
            entry = session.get(model, key)
            session.delete(entry)
            session.commit()
            settings.instance._update_cloud_sqlite_file()
            logger.success(
                f"Deleted {colors.yellow(f'row {key}')} in"
                f" {colors.blue(f'table {name}')}."
            )
            if name == "dobject":
                track_usage(entry.id, "delete")
                storage_key = storage_key_from_dobject(entry)
                delete_file(storage_key)
                logger.success(f"Deleted {colors.yellow(storage_key)} from storage.")

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
