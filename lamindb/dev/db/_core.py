from typing import Tuple

import sqlalchemy as sa
import sqlmodel as sqm
from lndb_setup import settings

from ...schema._table import Table


def session() -> sqm.Session:
    """Get connection session to DB engine.

    Returns a `sqlmodel.Session` object.
    """
    return sqm.Session(settings.instance.db_engine())


def get_foreign_keys(
    table_name: str, inspector=None, referred: Tuple[str, str] = None
) -> dict:
    """Return foreign keys of a table.

    Returns {constrained_column: (referred_table, referred_column)}
    """
    if inspector is None:
        inspector = sa.inspect(settings.instance.db_engine())

    keys = {}
    results = inspector.get_foreign_keys(table_name)
    if len(results) > 0:
        for result in results:
            referred_table = result["referred_table"]
            for i, j in zip(result["constrained_columns"], result["referred_columns"]):
                keys[i] = (referred_table, j)
    if referred is not None:
        keys = {k: v for k, v in results.items() if v == referred}
    return keys


def get_link_tables(inspector=None):
    """Link tables."""
    if inspector is None:
        inspector = sa.inspect(settings.instance.db_engine())
    link_tables = []
    for name in inspector.get_table_names():
        pks = inspector.get_pk_constraint(name)["constrained_columns"]
        columns = [i["name"] for i in inspector.get_columns(name)]
        if pks == columns and len(inspector.get_foreign_keys(name)) > 0:
            link_tables.append(name)

    return link_tables


def get_link_table(table1, table2):
    link_tables = get_link_tables()
    pks = Table.get_pks(table1) + Table.get_pks(table2)
    for table in link_tables:
        if set(pks) == set(Table.get_pks(table)):
            return table
    return None
