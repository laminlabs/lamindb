import sqlalchemy as sa
import sqlmodel as sqm
from lndb_setup import settings

from ...schema._table import Table


def session() -> sqm.Session:
    """Get connection session to DB engine.

    Returns a `sqlmodel.Session` object.
    """
    return sqm.Session(settings.instance.db_engine())


def get_foreign_keys(table_name: str, inspector=None, referred: tuple[str, str] = None):
    if inspector is None:
        inspector = sa.inspect(settings.instance.db_engine())
    result = {
        column["constrained_columns"][0]: (
            column["referred_table"],
            column["referred_columns"][0],
        )
        for column in inspector.get_foreign_keys(table_name)
    }
    if referred is not None:
        result = {k: v for k, v in result.items() if v == referred}
    return result


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
