import erdiagram
import sqlalchemy as sql
from lndb_setup import settings
from sqlalchemy import Column, Table


def draw(view=True):
    """Make a diagram of entity relationships."""
    metadata = get_db_metadata()
    graph = erdiagram.create_schema_graph(
        metadata=metadata,
        show_datatypes=False,
        show_indexes=False,  # ditto for indexes
        rankdir="TB",
        concentrate=False,  # Don't try to join the relation lines together
    )
    if view:
        erdiagram.view(graph)
    else:
        return graph


def list_entities():
    """Return all entities."""
    metadata = get_db_metadata()
    table_names = [table.name for table in metadata.sorted_tables]
    return table_names


def get_db_metadata():
    engine = settings.instance.db_engine()
    metadata = sql.MetaData(bind=engine)
    metadata.reflect()
    return metadata


def get_db_metadata_as_dict():
    metadata = get_db_metadata()
    return {
        table_name: get_table_metadata_as_dict(table)
        for table_name, table in metadata.tables.items()
    }


def get_table_metadata_as_dict(table: Table):
    return {
        "primary_keys": table.primary_key.columns.keys(),
        "foreign_keys": [
            (fk.column.table.key, fk.column.key) for fk in table.foreign_keys
        ],
        "columns": {
            column_name: get_column_metadata_as_dict(column)
            for column_name, column in table.columns.items()
        },
    }


def get_column_metadata_as_dict(column: Column):
    return {
        "key": column.key,
        "type": column.type,
        "primary_key": column.primary_key,
        "foreign_keys": column.foreign_keys,
        "nullable": column.nullable,
        "default": column.default,
    }
