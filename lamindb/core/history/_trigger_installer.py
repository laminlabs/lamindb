import enum
from abc import ABC, abstractmethod

from django.apps import apps
from django.db.backends.base.base import BaseDatabaseWrapper
from lamin_utils import logger


class HistoryEventTypes(enum.Enum):
    INSERT = 0
    UPDATE = 1
    DELETE = 2


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

    def update_history_triggers(self):
        tables = self._get_db_tables()

        tables_with_installed_triggers = self.get_tables_with_installed_triggers()

        tables_missing_triggers = tables.difference(tables_with_installed_triggers)

        for table in tables_missing_triggers:
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
    def __init__(self, connection: BaseDatabaseWrapper):
        super().__init__(connection=connection)

    def install_triggers(self, table: str):
        cursor = self.connection.cursor()

        # TODO: should we be extracting the primary key from info schema as well?
        # The assumption that all tables are keyed off 'id' seems brittle.

        # TODO: migration_history_id is only capturing the latest migration from
        # the `lamindb` app, and isn't taking bionty, etc into account.

        # TODO: populate update and delete triggers once the insert trigger works

        cursor.execute(f"""
CREATE OR REPLACE FUNCTION lamindb_history_{table}_insert_fn()
    RETURNS TRIGGER
    LANGUAGE PLPGSQL
AS $$
DECLARE latest_migration int := (SELECT
        CAST(SUBSTRING(name FROM '^(\\d+)') AS INTEGER)
    FROM django_migrations
    WHERE app = 'lamindb'
    ORDER BY CAST(SUBSTRING(name FROM '^(\\d+)') AS INTEGER) DESC
    LIMIT 1);
BEGIN
    INSERT INTO lamindb_history
    (id, migration_history_id, table_name, record_id, record_data, event_type, created_at)
    VALUES
    (gen_random_uuid(), latest_migration, '{table}', NEW.id, row_to_json(NEW), {HistoryEventTypes.INSERT.value}, now());
    RETURN NEW;
END;
$$
""")  # noqa: S608

        print("Trigger function installed")

        cursor.execute(f"""
CREATE OR REPLACE TRIGGER lamindb_history_{table}_insert_trigger
    AFTER INSERT
    ON {table}
    FOR EACH ROW
    EXECUTE PROCEDURE lamindb_history_{table}_insert_fn();
                       """)

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
