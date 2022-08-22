from lndb_setup import settings
from sqlmodel import Session

from ..schema._schema import alltables


def _create_delete_func(name: str, schema_module):
    def query_func(cls, id, **kwargs):
        with Session(settings.instance.db_engine()) as session:
            entry = session.get(schema_module, id)
            for k, v in kwargs.items():
                if isinstance(k, tuple):
                    k[1] = int(k[1])
                entry.__setattr__(k, v)
            session.delete(entry)
            session.commit()

    query_func.__name__ = name
    return query_func


class delete:
    """Delete an entry based on its primary identifier.

    Example:
    >>> delete.{entity}(id=1)
    """

    pass


for name, schema_module in alltables.items():
    func = _create_delete_func(name=name, schema_module=schema_module)
    setattr(delete, name, classmethod(func))
