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


FOREIGN_KEYS_LIST_COLUMN_NAME = "_lamin_fks"

RESERVED_COLUMNS: tuple[str] = (FOREIGN_KEYS_LIST_COLUMN_NAME,)

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
            if t not in existing_tables
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

    constraint_name: str
    constraint_type: Literal["PRIMARY KEY", "FOREIGN KEY"]

    # These need to be a list to account for composite primary keys
    source_columns: list[str]
    target_columns: list[str]

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
    ) -> tuple[KeyConstraint, list[KeyConstraint]]:
        cursor.execute(
            """
SELECT
    tc.conname AS constraint_name,
    CASE tc.contype
        WHEN 'p' THEN 'PRIMARY KEY'
        WHEN 'f' THEN 'FOREIGN KEY'
    END AS constraint_type,
    a.attname AS source_column,
    CASE WHEN tc.contype = 'f' THEN af.attname ELSE NULL END AS target_column,
    CASE WHEN tc.contype = 'f' THEN tf.relname ELSE NULL END AS target_table
FROM
    pg_constraint tc
    JOIN pg_class t ON t.oid = tc.conrelid
    JOIN pg_namespace n ON n.oid = t.relnamespace
    JOIN pg_attribute a ON a.attrelid = tc.conrelid AND a.attnum = ANY(tc.conkey)
    LEFT JOIN pg_class tf ON tf.oid = tc.confrelid
    LEFT JOIN pg_attribute af ON af.attrelid = tc.confrelid AND
        af.attnum = tc.confkey[array_position(tc.conkey, a.attnum)]
WHERE
    t.relname = %s
    AND tc.contype IN ('p', 'f')
ORDER BY
    tc.conname,
    array_position(tc.conkey, a.attnum);
""",
            (table,),
        )

        keys = cursor.fetchall()

        primary_key_constraint = None

        foreign_key_constraints: dict[str, KeyConstraint] = {}

        for k in keys:
            (
                constraint_name,
                constraint_type,
                source_column,
                target_column,
                target_table,
            ) = k

            if constraint_type == "PRIMARY KEY":
                if primary_key_constraint is None:
                    primary_key_constraint = KeyConstraint(
                        constraint_name=constraint_name,
                        constraint_type=constraint_type,
                        source_columns=[],
                        target_columns=[],
                        target_table=target_table,
                    )

                primary_key_constraint.source_columns.append(source_column)
                primary_key_constraint.target_columns.append(target_column)
            elif constraint_type == "FOREIGN KEY":
                if constraint_name not in foreign_key_constraints:
                    foreign_key_constraints[constraint_name] = KeyConstraint(
                        constraint_name=constraint_name,
                        constraint_type="FOREIGN KEY",
                        source_columns=[],
                        target_columns=[],
                        target_table=target_table,
                    )

                foreign_key_constraints[constraint_name].source_columns.append(
                    source_column
                )
                foreign_key_constraints[constraint_name].target_columns.append(
                    target_column
                )
            else:
                raise ValueError(f"Unhandled constraint type '{constraint_type}'")

        if primary_key_constraint is None:
            raise ValueError(
                f"Expected table {table} to have a primary key, but found none"
            )

        return (primary_key_constraint, list(foreign_key_constraints.values()))

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

        function_name = f"lamindb_history_{table}_fn"

        primary_key, foreign_key_constraints = self._get_table_key_constraints(
            table, cursor
        )

        table_columns = self._get_column_names(table, cursor)

        for col in table_columns:
            if col in RESERVED_COLUMNS:
                raise ValueError(
                    f"Column '{col}' in table {table} has a name that is "
                    "reserved by LaminDB for history tracking"
                )

        non_key_columns = table_columns.difference(set(primary_key.source_columns))

        for foreign_key_constraint in foreign_key_constraints:
            non_key_columns.difference_update(
                set(foreign_key_constraint.source_columns)
            )

        # We'll construct the record's data object from a list of keys and values.
        json_object_parts: list[str] = []

        # Foreign-key constraints will require the creation of variables within the function that
        # extract uniquely-identifiable columns from the record to which the foreign key refers.
        # We'll store the declarations of those variables here.
        foreign_key_uid_variables: dict[str, str] = {}

        # Any columns that are not part of a key can be added to the record's data directly.
        for column in non_key_columns:
            json_object_parts.extend([f"'{column}'", f"NEW.{column}"])

        for foreign_key_constraint in foreign_key_constraints:
            # Each foreign key constraint must be updated to refer to a set of uniquely identifying
            # columns within the target table.

            target_table_columns = self._get_column_names(
                foreign_key_constraint.target_table, cursor
            )

            uid_columns = []

            if "uid" in target_table_columns:
                uid_columns = ["uid"]
            else:
                # We don't know how to extract a unique identifier for this table, so we need to panic.
                raise ValueError(
                    f"Table '{table}' has a foreign key '{foreign_key_constraint.source_columns}' "
                    f"to '{foreign_key_constraint.target_table}'"
                    "that doesn't have a 'uid' column, so it's not clear what we need to "
                    f"serialize when storing {table}'s history"
                )

            # Give foreign-key variables names '_lamin_fkey_0', '_lamin_fkey_1', etc. to avoid
            # collisions. We'll store the table to which they foreign-key refers in the variable itself.
            fkey_variable_name = f"_lamin_fkey_{len(foreign_key_uid_variables)}"

            select_clause_parts = []
            null_select_clause_parts = []

            for uid_column in uid_columns:
                key_str = f"'{uid_column}'"
                select_clause_parts.extend([key_str, uid_column])
                null_select_clause_parts.extend([key_str, "NULL"])

            # Store the table's ID in historytablestate as part of the foreign-key data.
            source_column_name_array = ",".join(
                f"'{c}'" for c in foreign_key_constraint.source_columns
            )

            # We'll be creating a JSONB array with three elements:
            #  * the ID of the target table in HistoryTableState
            #  * the columns in the table that comprise the foreign-key constraint
            #  * the uniquely-identifiable columns in the target table that identify the target record
            where_clause = " AND ".join(
                f"{target_col} = NEW.{source_col}"
                for source_col, target_col in zip(
                    foreign_key_constraint.source_columns,
                    foreign_key_constraint.target_columns,
                )
            )

            fkey_variable_declaration = f"""
DECLARE {fkey_variable_name} jsonb :=
(
    SELECT
    jsonb_build_array(
        (SELECT id FROM lamindb_historytablestate WHERE table_name='{foreign_key_constraint.target_table}'),
        jsonb_build_array({source_column_name_array}),
        coalesce(
            (
                SELECT jsonb_build_object({", ".join(select_clause_parts)})
                FROM {foreign_key_constraint.target_table}
                WHERE {where_clause}
            ),
            jsonb_build_object({", ".join(null_select_clause_parts)})
        )
    )
);
"""  # noqa: S608

            foreign_key_uid_variables[fkey_variable_name] = fkey_variable_declaration

        # Place a list of foreign key constraints at the "end" of the object under a reserved name
        json_object_parts.extend(
            [
                f"'{FOREIGN_KEYS_LIST_COLUMN_NAME}'",
                f"jsonb_build_array({','.join(foreign_key_uid_variables.keys())})",
            ]
        )

        create_trigger_function_command = f"""
CREATE OR REPLACE FUNCTION {function_name}()
    RETURNS TRIGGER
    LANGUAGE PLPGSQL
AS $$
DECLARE latest_migration_id smallint := (
    SELECT MAX(id)
    FROM lamindb_historymigrationstate
);
DECLARE table_id smallint := (
    SELECT id
    FROM lamindb_historytablestate
    WHERE table_name='{table}'
    LIMIT 1
);
DECLARE history_triggers_locked bool := (SELECT EXISTS(SELECT locked FROM lamindb_historylock LIMIT 1) AND (SELECT locked FROM lamindb_historylock LIMIT 1));
{" ".join(foreign_key_uid_variables.values())}
DECLARE row_data jsonb;
DECLARE event_type int2;
DECLARE record_uid varchar(20);

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
            (id, history_migration_state_id, table_id, record_uid, record_data, event_type, created_at)
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
