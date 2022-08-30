from lndb_setup import settings
from sqlmodel import Session

from .._logger import colors, logger
from ..dev import storage_key_from_dobject, track_usage
from ..dev.file import delete_file
from ..schema._schema import alltables


def _create_delete_func(name: str, schema_module):
    def delete_func(cls, key):
        with Session(settings.instance.db_engine()) as session:
            entry = session.get(schema_module, key)
            session.delete(entry)
            session.commit()
            settings.instance._update_cloud_sqlite_file()
            logger.success(
                f"Deleted {colors.yellow(f'row {key}')} in"
                f" {colors.blue(f'table {name}')}."
            )
            if name == "dobject":
                track_usage(entry.id, entry.v, "delete")
                storage_key = storage_key_from_dobject(entry)
                delete_file(storage_key)
                logger.success(f"Deleted {colors.yellow(storage_key)} from storage.")

    delete_func.__name__ = name
    return delete_func


class delete:
    """Delete data by primary key.

    Deletes the row identified by the key in the table that represents the entity.
    For `dobject`, this also deletes the data object in storage.

    Example:

    >>> delete.dobject(key=("MHnsI0TWWNu3VBsedg6gS", "1"))
    """

    pass


for name, schema_module in alltables.items():
    func = _create_delete_func(name=name, schema_module=schema_module)
    setattr(delete, name, classmethod(func))
