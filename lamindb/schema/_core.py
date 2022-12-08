from typing import Union

import erdiagram
import sqlalchemy as sql
from lndb_setup import settings
from sqlalchemy import Column, Table


def view(save=False):
    """Make a diagram of entity relationships."""
    metadata = get_db_metadata()
    graph = erdiagram.create_schema_graph(
        metadata=metadata,
        show_datatypes=False,
        show_indexes=False,
        rankdir="TB",
        concentrate=True,
    )
    if not save:
        erdiagram.view(graph)
    else:
        return graph


def list_tables():
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
        "key": get_db_name(),
        "tables": {
            table_name: get_table_metadata_as_dict(table)
            for table_name, table in metadata.tables.items()
        },
    }


def get_table_metadata_as_dict(table: Table):
    return {
        "key": table.key,
        "primary_keys": table.primary_key.columns.keys(),
        "foreign_keys": get_foreign_keys_as_tuples(table),
        "columns": {
            column_name: get_column_metadata_as_dict(column)
            for column_name, column in table.columns.items()
        },
    }


def get_column_metadata_as_dict(column: Column):
    return {
        "key": column.key,
        "type": str(column.type),
        "primary_key": column.primary_key,
        "foreign_keys": get_foreign_keys_as_tuples(column),
        "nullable": column.nullable,
        "default": column.default,
    }


def get_foreign_keys_as_tuples(object: Union[Table, Column]):
    return [(fk.column.table.key, fk.column.key) for fk in object.foreign_keys]


def get_db_name() -> str:
    engine = settings.instance.db_engine()
    return engine.url.database


def get_table_object(table_name: str):
    metadata = get_db_metadata()
    return metadata.tables[table_name]
