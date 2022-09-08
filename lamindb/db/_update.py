from lndb_setup import settings
from sqlmodel import Session

from .._logger import colors, logger
from ..dev import track_usage
from ..schema._schema import alltables


def _create_update_func(name: str, schema_module):
    def update_func(cls, key, **kwargs):
        with Session(settings.instance.db_engine()) as session:
            entry = session.get(schema_module, key)
            for k, v in kwargs.items():
                if isinstance(k, tuple):
                    k[1] = int(k[1])
                entry.__setattr__(k, v)
            session.add(entry)
            session.commit()
            session.refresh(entry)
            logger.success(
                f"Updated entry {colors.green(f'{key}')} in {colors.blue(f'{name}')}."
            )
            if name == "dobject":
                track_usage(entry.id, entry.v, "update")

        settings.instance._update_cloud_sqlite_file()

    update_func.__name__ = name
    return update_func


class update:
    """Update an entry based on its primary identifier.

    Example:
    >>> update.{entity}(id=1, name='new_experiment')
    """

    pass


for name, schema_module in alltables.items():
    func = _create_update_func(name=name, schema_module=schema_module)
    setattr(update, name, classmethod(func))
