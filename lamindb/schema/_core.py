from typing import Union

import erdiagram
import sqlalchemy as sql
from lndb_setup._settings import SettingManager
from sqlalchemy import Column, Table


class Schema:
    def __init__(self, settings_manager: SettingManager) -> None:
        self.settings_manager = settings_manager

    def get_db_name(self) -> str:
        engine = self.settings_manager.instance.db_engine()
        return engine.url.database

    def get_db_metadata(self):
        engine = self.settings_manager.instance.db_engine()
        metadata = sql.MetaData(bind=engine)
        metadata.reflect()
        return metadata

    def draw(self, view=True):
        """Make a diagram of entity relationships."""
        metadata = self.get_db_metadata()
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

    def list_entities(self):
        """Return all entities."""
        metadata = self.get_db_metadata()
        table_names = [table.name for table in metadata.sorted_tables]
        return table_names

    def get_db_metadata_as_dict(self):
        metadata = self.get_db_metadata()
        return {
            "key": self.get_db_name(),
            "tables": {
                table_name: self.get_table_metadata_as_dict(table)
                for table_name, table in metadata.tables.items()
            },
        }

    def get_table_metadata_as_dict(self, table: Table):
        return {
            "key": table.key,
            "primary_keys": table.primary_key.columns.keys(),
            "foreign_keys": self.get_foreign_keys_as_tuples(table),
            "columns": {
                column_name: self.get_column_metadata_as_dict(column)
                for column_name, column in table.columns.items()
            },
        }

    def get_column_metadata_as_dict(self, column: Column):
        return {
            "key": column.key,
            "type": str(column.type),
            "primary_key": column.primary_key,
            "foreign_keys": self.get_foreign_keys_as_tuples(column),
            "nullable": column.nullable,
            "default": column.default,
        }

    def get_foreign_keys_as_tuples(self, object: Union[Table, Column]):
        return [(fk.column.table.key, fk.column.key) for fk in object.foreign_keys]

    def get_table_object(self, table_name: str):
        metadata = self.get_db_metadata()
        return metadata.tables[table_name]
