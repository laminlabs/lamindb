from typing import Union

import sqlmodel as sqm
from lndb_setup import settings

from ._select import SelectStmt, select


def get(
    table: sqm.SQLModel, primary_key=None, **fields
) -> Union[None, sqm.SQLModel, SelectStmt]:
    """Retrieve an entry from a table by primary key.

    This is a convenience function for retrieving rows from tables by fields of
    the table.

    It behaves like SQLAlchemy's/SQLModel's `session.get` for retrieving a row
    by primary key.

    In addition, it provides a shortcut for calling
    `select(table).where(table.col1 == val1, table.col2 == val2)` via
    `get(table, col1=val1, col2=val2)`.

    Args:
        table: Entity.
        primary_key: Primary key.
        fields: Fields and values passed as keyword arguments.
    """
    if primary_key is not None:
        if len(fields) > 0:
            raise RuntimeError("Either pass key or fields.")
        with sqm.Session(settings.instance.db_engine()) as session:
            return session.get(table, primary_key)
    elif len(fields) > 0:
        conditions = [getattr(table, k) == v for k, v in fields.items()]
        return select(table).where(*conditions)
    else:
        raise RuntimeError("Either pass primary_key or **fields.")
