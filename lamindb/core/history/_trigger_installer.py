import enum
import re
from abc import ABC, abstractmethod

from django.apps import apps
from django.db.backends.base.base import BaseDatabaseWrapper
from django.db.backends.utils import CursorWrapper
from lamin_utils import logger

from lamindb.models.history import HistoryLock


class HistoryEventTypes(enum.Enum):
    INSERT = 0
    UPDATE = 1
    DELETE = 2


# Certain tables must be excluded from history triggers to avoid
# creating infinite loops of triggering.
EXCLUDED_TABLES = ["lamindb_history", "lamindb_historylock"]


class HistoryRecordingTriggerInstaller(ABC):
    """Installs triggers that record history into the database's tables."""

    def __init__(self, connection: BaseDatabaseWrapper):
        self.connection = connection

    @abstractmethod
    def install_triggers(self, table: str):
        raise NotImplementedError()

    @abstractmethod
    def get_tables_with_installed_triggers(self) -> set[str]:
        raise NotImplementedError()

    def update_history_triggers(self, update_all: bool = False):
        # Ensure that the history lock exists
        HistoryLock.load()

        tables = self._get_db_tables()

        tables_with_installed_triggers = self.get_tables_with_installed_triggers()
        tables_missing_triggers = tables.difference(tables_with_installed_triggers)

        if update_all:
            tables_to_update = tables
        else:
            tables_to_update = tables_missing_triggers

        for table in tables_to_update:
            if table not in EXCLUDED_TABLES:
                logger.info(f"Installing history recording triggers for {table}")
                self.install_triggers(table)

    def _get_db_tables(self) -> set[str]:
        """Get the schema and table names for all tables in the database for which we might want to track history."""
        all_models = apps.get_models()

        tables = set()

        for model in all_models:
            # Record the model's underlying table, as well as
            # the underlying table for all of its many-to-many
            # relationships
            tables.add(model._meta.db_table)

            for field in model._meta.many_to_many:
                m2m_field = field
                through_model = m2m_field.remote_field.through

                tables.add(through_model._meta.db_table)

        return {t for t in tables if not t.startswith("django_")}


class PostgresHistoryRecordingTriggerInstaller(HistoryRecordingTriggerInstaller):
    # Since we're creating triggers and functions based on table names,
    # it's probably a good idea to check that those table names are valid
    # to mitigate the potential for SQL injection.
    VALID_TABLE_NAME_REGEX = re.compile("^[a-z_][a-z0-9_$]*$")

    def __init__(self, connection: BaseDatabaseWrapper):
        super().__init__(connection=connection)

    def _install_trigger(
        self,
        table: str,
        history_event_type: HistoryEventTypes,
        function_body_sql: str,
        cursor: CursorWrapper,
    ):
        function_name = f"lamindb_history_{table}_{history_event_type.name.lower()}_fn"
        trigger_name = f"lamindb_history_{table}_{history_event_type.name.lower()}_tr"

        cursor.execute(
            f"""
CREATE OR REPLACE FUNCTION {function_name}()
    RETURNS TRIGGER
    LANGUAGE PLPGSQL
AS $$
DECLARE latest_migration int := (SELECT
        CAST(SUBSTRING(name FROM '^(\\d+)') AS INTEGER)
    FROM django_migrations
    WHERE app = 'lamindb'
    ORDER BY CAST(SUBSTRING(name FROM '^(\\d+)') AS INTEGER) DESC
    LIMIT 1);
DECLARE history_triggers_locked bool := (SELECT EXISTS(SELECT locked FROM lamindb_historylock LIMIT 1) AND (SELECT locked FROM lamindb_historylock LIMIT 1));
BEGIN
    IF NOT history_triggers_locked THEN
        {function_body_sql}
    END IF;
    RETURN NEW;
END;
$$
""",  # noqa: S608
            (table, history_event_type.value),  # noqa: S608
        )

        cursor.execute(f"""
        CREATE OR REPLACE TRIGGER {trigger_name}
        AFTER {history_event_type.name}
        ON {table}
        FOR EACH ROW EXECUTE PROCEDURE {function_name}();
        """)

    def _install_insert_trigger(self, table: str, cursor: CursorWrapper) -> str:
        return self._install_trigger(
            table=table,
            history_event_type=HistoryEventTypes.INSERT,
            function_body_sql="""
        INSERT INTO lamindb_history
        (id, migration_history_id, table_name, record_id, record_data, event_type, created_at)
        VALUES
        (gen_random_uuid(), latest_migration, %s, NEW.id, row_to_json(NEW), %s, now());
        """,
            cursor=cursor,
        )

    def _install_update_trigger(self, table: str, cursor: CursorWrapper):
        return self._install_trigger(
            table=table,
            history_event_type=HistoryEventTypes.UPDATE,
            function_body_sql="""
        INSERT INTO lamindb_history
        (id, migration_history_id, table_name, record_id, record_data, event_type, created_at)
        VALUES
        (gen_random_uuid(), latest_migration, %s, NEW.id, row_to_json(NEW), %s, now());
        """,
            cursor=cursor,
        )

    def _install_delete_trigger(self, table: str, cursor: CursorWrapper):
        return self._install_trigger(
            table=table,
            history_event_type=HistoryEventTypes.DELETE,
            function_body_sql="""
        INSERT INTO lamindb_history
        (id, migration_history_id, table_name, record_id, record_data, event_type, created_at)
        VALUES
        (gen_random_uuid(), latest_migration, %s, OLD.id, NULL, %s, now());
        """,
            cursor=cursor,
        )

    def install_triggers(self, table: str):
        if not self.VALID_TABLE_NAME_REGEX.match(table):
            raise ValueError(
                f"Table name '{table}' doesn't look like a valid PostgreSQL table name"
            )

        cursor = self.connection.cursor()

        # TODO: should we be extracting the primary key from info schema as well?
        # The assumption that all tables are keyed off 'id' seems brittle.

        # TODO: migration_history_id is only capturing the latest migration from
        # the `lamindb` app, and isn't taking bionty, etc into account.

        self._install_insert_trigger(table, cursor)
        self._install_update_trigger(table, cursor)
        self._install_delete_trigger(table, cursor)

        print("Trigger functions installed")

    def get_tables_with_installed_triggers(self) -> set[str]:
        cursor = self.connection.cursor()

        cursor.execute("""
SELECT
    DISTINCT event_object_table
FROM information_schema.triggers
WHERE trigger_name like 'lamindb_history_%';
""")
        tables = set()

        rows = cursor.fetchall()

        for row in rows:
            table_name = row[0]
            tables.add(table_name)

        return tables


class SQLiteHistoryRecordingTriggerInstaller(HistoryRecordingTriggerInstaller):
    def __init__(self, connection: BaseDatabaseWrapper):
        super().__init__(connection=connection)

    def get_tables_with_installed_triggers(self) -> set[str]:
        logger.warning(
            "No triggers are installed since we're running against a SQLite backend"
        )
        return set()

    def install_triggers(self, table: str):
        logger.warning(
            f"Skipping trigger installation for table {'.'.join(table)} "
            "because we're running against a SQLite backend"
        )
        return


def create_history_recording_trigger_installer(
    connection: BaseDatabaseWrapper,
) -> HistoryRecordingTriggerInstaller:
    if connection.vendor == "postgresql":
        return PostgresHistoryRecordingTriggerInstaller(connection)
    elif connection.vendor == "sqlite":
        return SQLiteHistoryRecordingTriggerInstaller(connection)
    else:
        raise ValueError(f"History is not supported for vendor '{connection.vendor}'")
