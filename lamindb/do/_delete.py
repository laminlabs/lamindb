from lndb_setup import settings
from sqlmodel import Session

from .._logger import colors, logger
from ..dev import track_usage
from ..schema._schema import alltables


def _create_delete_func(name: str, schema_module):
    def delete_func(cls, id, **kwargs):
        with Session(settings.instance.db_engine()) as session:
            entry = session.get(schema_module, id)
            for k, v in kwargs.items():
                if isinstance(k, tuple):
                    k[1] = int(k[1])
                entry.__setattr__(k, v)
            session.delete(entry)
            session.commit()
            logger.success(
                f"Deleted {colors.yellow(f'entry {entry.id}')} in"
                f" {colors.blue(f'table {name}')}!"
            )
            if name == "dobject":
                track_usage(entry.id, entry.v, "delete")

        settings.instance._update_cloud_sqlite_file()

    delete_func.__name__ = name
    return delete_func


class delete:
    """Delete an entry based on its primary identifier.

    Example:
    >>> delete.{entity}(id=1)
    """

    pass


for name, schema_module in alltables.items():
    func = _create_delete_func(name=name, schema_module=schema_module)
    setattr(delete, name, classmethod(func))
