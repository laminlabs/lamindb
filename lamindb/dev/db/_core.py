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


def check_if_link_table(table_name: str):
    """Check if a table is a link table.

    We define link tables there is overlap between primary and foreign keys
    """
    pks = Table.get_pks(table_name)
    fks = get_foreign_keys(table_name).keys()
    intersect = set(pks).intersection(fks)
    if intersect:
        return intersect


def get_link_tables(inspector=None):
    """Get all link tables."""
    if inspector is None:
        inspector = sa.inspect(settings.instance.db_engine())
    link_tables = []
    for name in inspector.get_table_names():
        if check_if_link_table(name):
            link_tables.append(name)

    return link_tables


def get_link_table(table1: str, table2: str):
    link_tables = get_link_tables()
    pks = [f"{table1}_{i}" for i in Table.get_pks(table1)] + [
        f"{table2}_{i}" for i in Table.get_pks(table2)
    ]

    for table in link_tables:
        if set(pks) == set(Table.get_pks(table)):
            return table
    return None
