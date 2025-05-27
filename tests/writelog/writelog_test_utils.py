from django.db.backends.utils import CursorWrapper
from lamindb.core.writelog._db_metadata_wrapper import PostgresDatabaseMetadataWrapper
from lamindb.core.writelog._types import UIDColumns
from typing_extensions import override


class FakeMetadataWrapper(PostgresDatabaseMetadataWrapper):
    """A fake DB metadata wrapper that allows us to control which database tables the installer will see and target."""

    def __init__(self):
        super().__init__()
        self._tables_with_triggers = set()
        self._db_tables: set[str] = set()
        self._many_to_many_tables: set[str] = set()
        self._uid_columns: dict[str, UIDColumns] = {}

    @override
    def get_tables_with_installed_triggers(self, cursor: CursorWrapper) -> set[str]:
        return self._tables_with_triggers

    def set_tables_with_installed_triggers(self, tables: set[str]):
        self._tables_with_triggers = tables

    @override
    def get_db_tables(self) -> set[str]:
        return self._db_tables

    def set_db_tables(self, tables: set[str]):
        self._db_tables = tables

    @override
    def get_many_to_many_db_tables(self) -> set[str]:
        return self._many_to_many_tables

    def set_many_to_many_db_tables(self, tables: set[str]):
        self._many_to_many_tables = tables

    @override
    def get_uid_columns(self, table: str, cursor: CursorWrapper) -> UIDColumns:
        if table in self._uid_columns:
            return self._uid_columns[table]
        else:
            return super().get_uid_columns(table, cursor)

    def set_uid_columns(self, table: str, uid_columns: UIDColumns):
        self._uid_columns[table] = uid_columns
