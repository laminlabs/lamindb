import enum
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Literal

from django.apps import apps
from django.db import transaction
from django.db.backends.base.base import BaseDatabaseWrapper
from django.db.backends.utils import CursorWrapper
from lamin_utils import logger

from lamindb.models.history import HistoryLock, HistoryTableState


class HistoryEventTypes(enum.Enum):
    INSERT = 0
    UPDATE = 1
    DELETE = 2


# Certain tables must be excluded from history triggers to avoid
# creating infinite loops of triggering. Others we're excluding
# because their state is managed by LaminDB's underlying Django
# machinery.
EXCLUDED_TABLES = [
    "lamindb_history",
    "lamindb_historylock",
    "lamindb_historytablestate",
    "lamindb_historymigrationstate",
    "django_content_type",
    "django_migrations",
    # Skipping these temporarily
    "lamindb_schemaparam",
    "lamindb_artifactparamvalue",
    "lamindb_param",
    "lamindb_paramvalue",
    "lamindb_runparamvalue",
    "lamindb_rundata",
    "lamindb_artifactfeaturevalue",
    "lamindb_flextabledata",
]


class HistoryRecordingTriggerInstaller(ABC):
    """Installs triggers that record history into the database's tables."""

    def __init__(self, connection: BaseDatabaseWrapper):
        self.connection = connection

    @abstractmethod
    def install_triggers(self, table: str, cursor: CursorWrapper):
        raise NotImplementedError()

    @abstractmethod
    def get_tables_with_installed_triggers(self, cursor: CursorWrapper) -> set[str]:
        raise NotImplementedError()

    def _update_history_table_state(self, tables: set[str]):
        existing_tables = set(
            HistoryTableState.objects.filter(table_name__in=tables).values_list(
                "table_name", flat=True
            )
        )

        table_states_to_add = [
            HistoryTableState(table_name=t, backfilled=False)
            for t in tables
            if tables not in existing_tables
        ]

        HistoryTableState.objects.bulk_create(table_states_to_add)

    def _update_migration_state(self, cursor: CursorWrapper):
        cursor.execute("""
WITH
-- Get current migration state
current_state AS (
    SELECT jsonb_agg(t) AS state FROM (
        SELECT
            MAX(CAST(SUBSTRING(name FROM '^([0-9]+)_') AS INTEGER)) as migration_id,
            app
        FROM django_migrations
        GROUP BY app
    ) t
),
-- Get latest stored state
latest_state AS (
    SELECT migration_history_id AS state
    FROM lamindb_historymigrationstate
    ORDER BY id DESC
    LIMIT 1
)
INSERT INTO lamindb_historymigrationstate (migration_history_id)
SELECT current_state.state
FROM current_state
WHERE NOT EXISTS (
    SELECT 1 FROM latest_state WHERE latest_state.state = current_state.state
)
AND current_state.state IS NOT NULL;
""")

    def update_history_triggers(self, update_all: bool = False):
        with transaction.atomic():
            # Ensure that the history lock exists
            HistoryLock.load()
            tables = self._get_db_tables()
            self._update_history_table_state(tables)

            cursor = self.connection.cursor()

            self._update_migration_state(cursor)

            tables_with_installed_triggers = self.get_tables_with_installed_triggers(
                cursor
            )
            tables_missing_triggers = tables.difference(tables_with_installed_triggers)

            if update_all:
                tables_to_update = tables
            else:
                tables_to_update = tables_missing_triggers

            for table in tables_to_update:
                if table not in EXCLUDED_TABLES:
                    logger.important(
                        f"Installing history recording triggers for {table}"
                    )
                    self.install_triggers(table, cursor)

    def _get_db_tables(self) -> set[str]:
        """Get the schema and table names for all tables in the database for which we might want to track history."""
        all_models = apps.get_models()

        tables: set[str] = set()

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


@dataclass
class KeyConstraint:
    """Simple encapsulation of one of a table's key constraints."""

    constraint_type: Literal["PRIMARY KEY", "FOREIGN KEY"]
    source_column: str
    target_table: str


class PostgresHistoryRecordingTriggerInstaller(HistoryRecordingTriggerInstaller):
    # Since we're creating triggers and functions based on table names,
    # it's probably a good idea to check that those table names are valid
    # to mitigate the potential for SQL injection.
    VALID_TABLE_NAME_REGEX = re.compile("^[a-z_][a-z0-9_$]*$")

    def __init__(self, connection: BaseDatabaseWrapper):
        super().__init__(connection=connection)

    def _get_table_key_constraints(
        self, table: str, cursor: CursorWrapper
    ) -> list[KeyConstraint]:
        cursor.execute(
            """
SELECT
    tc.constraint_type,
    kcu.column_name AS source_column,
    ccu.table_name AS target_table
FROM
    information_schema.table_constraints AS tc
    JOIN information_schema.key_column_usage AS kcu
      ON tc.constraint_name = kcu.constraint_name
      AND tc.table_schema = kcu.table_schema
    JOIN information_schema.constraint_column_usage AS ccu
      ON ccu.constraint_name = tc.constraint_name
      AND ccu.table_schema = tc.table_schema
WHERE
    tc.table_name = %s and constraint_type in ('PRIMARY KEY', 'FOREIGN KEY');
""",
            (table,),
        )

        keys = cursor.fetchall()

        return [
            KeyConstraint(constraint_type=k[0], source_column=k[1], target_table=k[2])
            for k in keys
        ]

    def _get_column_names(self, table: str, cursor: CursorWrapper) -> set[str]:
        cursor.execute(
            "SELECT column_name FROM information_schema.columns WHERE TABLE_NAME = %s ORDER BY ordinal_position",
            (table,),
        )

        return {r[0] for r in cursor.fetchall()}

    def install_triggers(self, table: str, cursor: CursorWrapper):
        if not self.VALID_TABLE_NAME_REGEX.match(table):
            raise ValueError(
                f"Table name '{table}' doesn't look like a valid PostgreSQL table name"
            )

        # TODO: migration_history_id is only capturing the latest migration from
        # the `lamindb` app, and isn't taking bionty, etc into account.

        function_name = f"lamindb_history_{table}_fn"

        key_constraints = self._get_table_key_constraints(table, cursor)

        primary_keys: set[str] = {
            c.source_column
            for c in key_constraints
            if c.constraint_type == "PRIMARY KEY"
        }
        foreign_key_constraints: dict[str, KeyConstraint] = {
            c.source_column: c
            for c in key_constraints
            if c.constraint_type == "FOREIGN KEY"
        }

        table_columns = self._get_column_names(table, cursor)

        json_object_parts: list[str] = []
        foreign_key_uid_variables = []

        for column in table_columns.difference(primary_keys):
            if column in foreign_key_constraints:
                key_constraint = foreign_key_constraints[column]

                if "uid" not in self._get_column_names(
                    key_constraint.target_table, cursor
                ):
                    raise ValueError(
                        f"Table '{table}' has a foreign key to '{key_constraint.target_table}' "
                        "that doesn't have a 'uid' column, so it's not clear what we need to "
                        f"serialize when storing {table}'s history"
                    )

                primary_key_lookup_clauses = " AND ".join(
                    f"{c} = NEW.{c}" for c in primary_keys
                )
                foreign_key_uid_var = f"fkey_{column}_uid"

                foreign_key_uid_variables.append(f"""
DECLARE {foreign_key_uid_var} varchar(8) :=
(
    SELECT uid FROM {key_constraint.target_table}
    WHERE {primary_key_lookup_clauses}
);
                    """)  # noqa: S608
                json_object_parts.extend([f"'{column}._uid'", foreign_key_uid_var])
            else:
                json_object_parts.extend([f"'{column}'", f"NEW.{column}"])

        create_trigger_function_command = f"""
CREATE OR REPLACE FUNCTION {function_name}()
    RETURNS TRIGGER
    LANGUAGE PLPGSQL
AS $$
DECLARE latest_migration_id smallint := (
    SELECT MAX(migration_history_id)
    FROM lamindb_historymigrationstate
);
DECLARE table_id smallint := (
    SELECT id
    FROM lamindb_historytablestate
    WHERE table_name='{table}'
    LIMIT 1
);
DECLARE history_triggers_locked bool := (SELECT EXISTS(SELECT locked FROM lamindb_historylock LIMIT 1) AND (SELECT locked FROM lamindb_historylock LIMIT 1));
{" ".join(foreign_key_uid_variables)}
DECLARE row_data jsonb;
DECLARE event_type int2;
DECLARE record_uid varchar;

BEGIN
    IF NOT history_triggers_locked THEN

        IF (TG_OP = 'DELETE') THEN
            row_data := NULL;
            record_uid := {"OLD.uid" if "uid" in table_columns else "NULL"};
            event_type := {HistoryEventTypes.DELETE.value};
        ELSE
            row_data := jsonb_build_object({", ".join(json_object_parts)});
            record_uid := {"NEW.uid" if "uid" in table_columns else "NULL"};

            IF (TG_OP = 'INSERT') THEN
                event_type := {HistoryEventTypes.INSERT.value};
            ELSIF (TG_OP = 'UPDATE') THEN
                event_type := {HistoryEventTypes.UPDATE.value};
            END IF;
        END IF;

        INSERT INTO lamindb_history
            (id, migration_history_id, table_id, record_uid, record_data, event_type, created_at)
        VALUES
            (gen_random_uuid(), latest_migration_id, table_id, record_uid, row_data, event_type, now());
    END IF;
    RETURN NEW;
END;
$$
"""  # noqa: S608
        logger.debug(create_trigger_function_command)
        cursor.execute(create_trigger_function_command)

        for history_event_type in HistoryEventTypes:
            trigger_name = (
                f"lamindb_history_{table}_{history_event_type.name.lower()}_tr"
            )

            cursor.execute(f"""
            CREATE OR REPLACE TRIGGER {trigger_name}
            AFTER {history_event_type.name}
            ON {table}
            FOR EACH ROW EXECUTE PROCEDURE {function_name}();
            """)

    def get_tables_with_installed_triggers(self, cursor: CursorWrapper) -> set[str]:
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

    def get_tables_with_installed_triggers(self, cursor: CursorWrapper) -> set[str]:
        logger.warning(
            "No triggers are installed since we're running against a SQLite backend"
        )
        return set()

    def install_triggers(self, table: str, cursor: CursorWrapper):
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
