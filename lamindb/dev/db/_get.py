from typing import Union, overload  # noqa

import sqlmodel as sqm
from lndb_setup import settings


def get(table: sqm.SQLModel, primary_key) -> Union[None, sqm.SQLModel]:
    """Retrieve an entry from a table by primary key or fields.

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
    """
    with sqm.Session(settings.instance.db_engine()) as session:
        return session.get(table, primary_key)
