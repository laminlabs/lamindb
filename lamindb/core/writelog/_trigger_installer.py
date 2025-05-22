import enum
import re
from abc import ABC, abstractmethod
from typing import Any

from django.db import models, transaction
from django.db.backends.base.base import BaseDatabaseWrapper
from django.db.backends.utils import CursorWrapper
from lamin_utils import logger
from typing_extensions import override

from lamindb.models.writelog import (
    DEFAULT_BRANCH_CODE,
    DEFAULT_CREATED_BY_UID,
    DEFAULT_RUN_UID,
    WriteLogLock,
    WriteLogMigrationState,
    WriteLogTableState,
)

from ._db_metadata_wrapper import (
    DatabaseMetadataWrapper,
    PostgresDatabaseMetadataWrapper,
)
from ._types import KeyConstraint, UIDColumns


class WriteLogEventTypes(enum.Enum):
    INSERT = 0
    UPDATE = 1
    DELETE = 2


class DjangoMigration(models.Model):
    """This model class allows us to access the migrations table using normal Django syntax."""

    app = models.CharField(max_length=255)
    name = models.CharField(max_length=255)
    applied = models.DateTimeField()

    class Meta:
        managed = False  # Tell Django not to manage this table
        db_table = "django_migrations"  # Specify the actual table name


FOREIGN_KEYS_LIST_COLUMN_NAME = "_lamin_fks"

RESERVED_COLUMNS: tuple[str] = (FOREIGN_KEYS_LIST_COLUMN_NAME,)

# Certain tables must be excluded from write log triggers to avoid
# creating infinite loops of triggering. Others we're excluding
# because their state is managed by LaminDB's underlying Django
# machinery.
EXCLUDED_TABLES = [
    "lamindb_writelog",
    "lamindb_writeloglock",
    "lamindb_writelogtablestate",
    "lamindb_writelogmigrationstate",
    "django_content_type",
    "django_migrations",
    # FIXME what to do with these?
    "lamindb_flextabledata",
    "lamindb_rundata",
]


class WriteLogRecordingTriggerInstaller(ABC):
    """Installs triggers that record write logs into the database's tables."""

    def __init__(
        self,
        connection: BaseDatabaseWrapper,
        db_metadata: DatabaseMetadataWrapper,
    ):
        self.connection = connection
        self.db_metadata = db_metadata

    @abstractmethod
    def install_triggers(self, table: str, cursor: CursorWrapper):
        raise NotImplementedError()

    def _update_write_log_table_state(self, tables: set[str]):
        existing_tables = set(
            WriteLogTableState.objects.filter(table_name__in=tables).values_list(
                "table_name", flat=True
            )
        )

        table_states_to_add = [
            WriteLogTableState(table_name=t, backfilled=False)
            for t in tables
            if t not in existing_tables
        ]

        WriteLogTableState.objects.bulk_create(table_states_to_add)

    def _update_migration_state(self, cursor: CursorWrapper):
        app_migrations = {}

        for app in DjangoMigration.objects.values_list("app", flat=True).distinct():
            migrations = DjangoMigration.objects.filter(app=app)
            max_migration_id = 0

            for migration in migrations:
                # Extract the number from the migration name
                match = re.match(r"^([0-9]+)_", migration.name)  # type: ignore
                if match:
                    migration_id = int(match.group(1))
                    max_migration_id = max(max_migration_id, migration_id)

            app_migrations[app] = max_migration_id

        current_state = [
            {"migration_id": mig_id, "app": app}
            for app, mig_id in sorted(app_migrations.items())
        ]

        try:
            latest_state = WriteLogMigrationState.objects.order_by("-id").first()
            latest_state_json = (
                latest_state.migration_state_id if latest_state else None
            )
        except WriteLogMigrationState.DoesNotExist:
            latest_state_json = None

        if current_state and current_state != latest_state_json:
            WriteLogMigrationState.objects.create(migration_state_id=current_state)

    def update_write_log_triggers(self, update_all: bool = False):
        with transaction.atomic():
            # Ensure that the write log lock exists
            WriteLogLock.load()
            tables = self.db_metadata.get_db_tables()
            self._update_write_log_table_state(tables)

            cursor = self.connection.cursor()

            self._update_migration_state(cursor)

            tables_with_installed_triggers = (
                self.db_metadata.get_tables_with_installed_triggers(cursor)
            )
            tables_missing_triggers = tables.difference(tables_with_installed_triggers)

            if update_all:
                tables_to_update = tables
            else:
                tables_to_update = tables_missing_triggers

            for table in tables_to_update:
                if table not in EXCLUDED_TABLES:
                    logger.important(
                        f"installing write log recording triggers for {table}"
                    )
                    self.install_triggers(table, cursor)


class PostgresTriggerBuilder:
    # Since we're creating triggers and functions based on table names,
    # it's probably a good idea to check that those names are valid
    # to mitigate the potential for SQL injection.
    VALID_TABLE_NAME_REGEX = re.compile("^[a-z_][a-z0-9_$]*$")
    VALID_VARIABLE_NAME_REGEX = re.compile("^[a-z_0-9]+$")

    def __init__(
        self, table: str, db_metadata: DatabaseMetadataWrapper, cursor: CursorWrapper
    ):
        self.table = table
        self.function_name = f"lamindb_writelog_{table}_fn"
        self.db_metadata = db_metadata
        self.cursor = cursor

        self._variables: dict[str, tuple[str, str | None]] = {}
        self._variable_order: list[str] = []

        self._record_uid: dict[bool, str] = {}

    def declare_variable(self, name: str, type: str, declaration: str | None = None):
        if not self.VALID_VARIABLE_NAME_REGEX.match(name):
            raise ValueError(f"'{name}' is not a valid Postgres variable name")

        if type not in ("int", "int2", "smallint", "bool", "jsonb", "varchar"):
            raise ValueError(
                f"Unknown variable type '{type}'. declare_variable()'s known type list is not "
                "comprehensive; if this is a valid type, add it to that function's type enum"
            )

        if name not in self._variables:
            self._variables[name] = (type, declaration)
            self._variable_order.append(name)

    def _build_variable_declaration(self, name: str) -> str:
        variable_type, declaration = self._variables[name]

        if declaration is None:
            return f"DECLARE {name} {variable_type};"
        else:
            return f"DECLARE {name} {variable_type} := ({declaration});"

    def _build_variable_declarations(self) -> str:
        """Build variable declarations in the order in which they were declared."""
        return "\n".join(
            self._build_variable_declaration(name) for name in self._variable_order
        )

    def _build_jsonb_array(self, items: list[str]) -> str:
        return f"jsonb_build_array({', '.join(items)})"

    def _build_jsonb_object(
        self, obj: dict[str | int, Any], escape_keys: bool = True
    ) -> str:
        parts = []

        for key, value in obj.items():
            if escape_keys and not isinstance(key, int):
                key = f"'{key}'"

            parts.append(key)
            parts.append(value)

        return f"jsonb_build_object({', '.join(str(x) for x in parts)})"

    def _build_record_uid(self, is_delete: bool) -> str:
        # Cache record UID computation
        if is_delete not in self._record_uid:
            self._record_uid[is_delete] = self._build_record_uid_inner(is_delete)

        return self._record_uid[is_delete]

    def _build_space_uid(self, is_delete: bool) -> str:
        columns = self.db_metadata.get_column_names(
            table=self.table, cursor=self.cursor
        )

        if "space_id" in columns:
            # All tables that are not many-to-many tables will have a space_id column
            if is_delete:
                table_name_in_trigger = "OLD"
            else:
                table_name_in_trigger = "NEW"

            return f"SELECT uid FROM lamindb_space WHERE id = {table_name_in_trigger}.space_id"  # noqa: S608
        else:
            # For the rest, we can set this to NULL.
            return "NULL"

    def _build_record_uid_inner(self, is_delete: bool) -> str:
        uid_columns_list: UIDColumns = self.db_metadata.get_uid_columns(
            table=self.table, cursor=self.cursor
        )

        if len(uid_columns_list) == 1:
            table_uid = uid_columns_list[0]

            if table_uid.source_table_name != self.table:
                raise ValueError(
                    f"Expected UID column for standard table {self.table} "
                    f"to refer to its own columns, but it refers to those of table {table_uid.source_table_name}"
                )

            # This table's UID is determined solely by its own columns
            if is_delete:
                table_name_in_trigger = "OLD"
            else:
                table_name_in_trigger = "NEW"

            return self._build_jsonb_array(
                [
                    f"{table_name_in_trigger}.{column}"
                    for column in table_uid.uid_columns
                ]
            )
        else:
            # This table's UID is the composition of the UIDs of the tables it's referencing.
            # Since it's possible to have multiple foreign-key references to the same table,
            # we'll store the UID as a list of three-tuples, each of which is
            # [foreign table ID, list of columns for the foreign key constraint, UID fields for the foreign-keyed table].
            record_uid: list[str] = []

            for table_uid in uid_columns_list:
                table_id_var = self.add_table_id_variable(table_uid.source_table_name)

                foreign_key_constraint = table_uid.key_constraint

                if foreign_key_constraint is None:
                    raise ValueError(
                        f"Expected UID for '{self.table}'s referent '{table_uid.source_table_name}' to "
                        "have an associated key constraint, but it didn't"
                    )

                fkey_uid_lookup_variable = self.add_foreign_key_uid_lookup_variable(
                    foreign_key_constraint=foreign_key_constraint, is_delete=is_delete
                )

                record_uid.append(
                    self._build_jsonb_array(
                        [
                            table_id_var,
                            self._build_jsonb_array(
                                [
                                    f"'{c}'"
                                    for c in foreign_key_constraint.source_columns
                                ]
                            ),
                            fkey_uid_lookup_variable,
                        ]
                    )
                )

            return self._build_jsonb_array(record_uid)

    def _build_record_data(self) -> str:
        table_columns = self.db_metadata.get_column_names(self.table, self.cursor)
        uid_columns_list = self.db_metadata.get_uid_columns(self.table, self.cursor)
        primary_key, foreign_key_constraints = (
            self.db_metadata.get_table_key_constraints(self.table, self.cursor)
        )

        # Make sure the table isn't using any reserved column names.
        for col in table_columns:
            if col in RESERVED_COLUMNS:
                raise ValueError(
                    f"Column '{col}' in table {self.table} has a name that is "
                    "reserved by LaminDB for write log tracking"
                )

        # We don't need to store the table's primary key columns in its data object,
        # since the object will be identified by its UID columns.
        non_key_columns = table_columns.difference(set(primary_key.source_columns))

        # We also don't need to store any foreign-key columns in its data object,
        # since the references those columns encode will be captured by reference to
        # their UID columns.
        for foreign_key_constraint in foreign_key_constraints:
            non_key_columns.difference_update(
                set(foreign_key_constraint.source_columns)
            )

        # Since we're recording the record's UID, we don't need to store the UID columns in the
        # data object as well.
        if len(uid_columns_list) == 1:
            table_uid = uid_columns_list[0]
            non_key_columns.difference_update(table_uid.uid_columns)

        record_data: dict[str | int, Any] = {}

        # Any columns that are not part of a key can be added to the record's data directly.
        for column in non_key_columns:
            record_data[column] = f"NEW.{column}"

        # Place a list of foreign key constraints at the "end" of the object under a reserved name
        fkey_source_columns_variables = [
            self.add_foreign_key_source_columns_variable(
                foreign_key_constraint=foreign_key_constraint, is_delete=False
            )
            for foreign_key_constraint in foreign_key_constraints
            # Don't record foreign-keys to space, since we store space_uid separately
            if not (
                foreign_key_constraint.target_table == "lamindb_space"
                and foreign_key_constraint.source_columns == ["space_id"]
            )
        ]

        record_data[FOREIGN_KEYS_LIST_COLUMN_NAME] = self._build_jsonb_array(
            fkey_source_columns_variables
        )

        return self._build_jsonb_object(record_data)

    def add_table_id_variable(self, table: str) -> str:
        """Adds a variable looking up a table's ID given its name.

        Will only add a variable for a given table name once.
        """
        table_id = WriteLogTableState.objects.get(table_name=table).id

        variable_name = f"_table_name_{table_id}"

        if variable_name not in self._variables:
            self._validate_table_name(table)

            self.declare_variable(
                variable_name,
                "smallint",
                f"SELECT id FROM lamindb_writelogtablestate WHERE table_name = '{table}' LIMIT 1",  # noqa: S608
            )

        return variable_name

    def add_foreign_key_source_columns_variable(
        self, foreign_key_constraint: KeyConstraint, is_delete: bool
    ) -> str:
        # Foreign-key constraints will require the creation of variables within the function that
        # extract uniquely-identifiable columns from the record to which the foreign key refers.
        # We'll store the declarations of those variables here.

        table_id = self.add_table_id_variable(foreign_key_constraint.target_table)

        uid_lookup_variable = self.add_foreign_key_uid_lookup_variable(
            foreign_key_constraint, is_delete=is_delete
        )

        # Give foreign-key variables names '_lamin_fkey_ref_0', '_lamin_fkey_ref_1', etc. to avoid
        # collisions. We'll store the table to which they foreign-key refers in the variable itself.
        variable_name = self._create_unique_variable_name("_lamin_fkey_ref")

        # We'll be creating a JSONB array with three elements:
        #  * the ID of the target table in WriteLogTableState
        #  * the columns in the table that comprise the foreign-key constraint
        #  * the uniquely-identifiable columns in the target table that identify the target record
        source_columns_lookup_array = self._build_jsonb_array(
            [
                table_id,
                self._build_jsonb_array(
                    [f"'{c}'" for c in foreign_key_constraint.source_columns]
                ),
                uid_lookup_variable,
            ]
        )

        self.declare_variable(
            name=variable_name,
            type="jsonb",
            declaration=f"SELECT {source_columns_lookup_array}",
        )

        return variable_name

    def _create_unique_variable_name(self, prefix: str) -> str:
        var_sequence_number = len(
            [v for v in self._variables.keys() if v.startswith(prefix)]
        )
        variable_name = f"{prefix}_{var_sequence_number}"

        return variable_name

    def add_foreign_key_uid_lookup_variable(
        self, foreign_key_constraint: KeyConstraint, is_delete: bool
    ) -> str:
        source_record = "OLD" if is_delete else "NEW"

        uid_column_list = self.db_metadata.get_uid_columns(
            table=foreign_key_constraint.target_table, cursor=self.cursor
        )

        if len(uid_column_list) > 1:
            raise ValueError(
                f"Table '{self.table}' has a foreign-key relationship "
                f"with '{foreign_key_constraint.target_table}', which appears to be a "
                "many-to-many table. This is not currently supported."
            )

        table_uid = uid_column_list[0]

        where_clause = " AND ".join(
            f"{target_col} = {source_record}.{source_col}"
            for source_col, target_col in zip(
                foreign_key_constraint.source_columns,
                foreign_key_constraint.target_columns,
            )
        )

        variable_name = self._create_unique_variable_name(
            prefix="_fkey_uid_d" if is_delete else "_fkey_uid"
        )

        self._validate_table_name(foreign_key_constraint.target_table)

        self.declare_variable(
            variable_name,
            "jsonb",
            f"""
coalesce(
    (
        SELECT {self._build_jsonb_object({c: c for c in table_uid.uid_columns})}
        FROM {foreign_key_constraint.target_table}
        WHERE {where_clause}
    ),
    {self._build_jsonb_object(dict.fromkeys(table_uid.uid_columns, "NULL"))}
)
""",  # noqa: S608
        )

        return variable_name

    def _validate_table_name(self, table_name: str):
        if not self.VALID_TABLE_NAME_REGEX.match(table_name):
            raise ValueError(
                f"Table name '{table_name}' doesn't look like a valid PostgreSQL table name"
            )

    def build(self) -> str:
        self._validate_table_name(self.table)

        self.declare_variable(
            name="write_log_triggers_locked",
            type="bool",
            declaration="SELECT EXISTS(SELECT locked FROM lamindb_writeloglock LIMIT 1) AND "
            "(SELECT locked FROM lamindb_writeloglock LIMIT 1)",
        )

        self.declare_variable(
            name="latest_migration_id",
            type="smallint",
            declaration="SELECT MAX(id) FROM lamindb_writelogmigrationstate",
        )

        self.declare_variable(name="space_uid", type="varchar")

        table_id_var = self.add_table_id_variable(table=self.table)
        self.declare_variable(name="record_data", type="jsonb")
        self.declare_variable(name="event_type", type="int2")
        self.declare_variable(name="record_uid", type="jsonb")

        # Build the function body so that any variables created as a side-effect of constructing
        # record data get populated.

        function_body = f"""
    IF NOT write_log_triggers_locked THEN
        IF (TG_OP = 'DELETE') THEN
            record_data := NULL;
            record_uid := {self._build_record_uid(is_delete=True)};
            event_type := {WriteLogEventTypes.DELETE.value};
            space_uid := ({self._build_space_uid(is_delete=True)});
        ELSE
            record_data := {self._build_record_data()};
            record_uid := {self._build_record_uid(is_delete=False)};
            space_uid := ({self._build_space_uid(is_delete=False)});

            IF (TG_OP = 'INSERT') THEN
                event_type := {WriteLogEventTypes.INSERT.value};
            ELSIF (TG_OP = 'UPDATE') THEN
                event_type := {WriteLogEventTypes.UPDATE.value};
            END IF;
        END IF;

        INSERT INTO lamindb_writelog
            (
                uid,
                migration_state_id,
                table_id,
                record_uid,
                space_uid,
                created_by_uid,
                branch_code,
                run_uid,
                record_data,
                event_type,
                created_at
            )
        VALUES
            (
                base62(18),
                latest_migration_id,
                {table_id_var},
                record_uid,
                space_uid,
                '{DEFAULT_CREATED_BY_UID}',
                {DEFAULT_BRANCH_CODE},
                '{DEFAULT_RUN_UID}',
                record_data,
                event_type,
                now()
            );
    END IF;
    RETURN NEW;
"""  # noqa: S608

        # Build the actual function declaration
        return f"""
CREATE OR REPLACE FUNCTION {self.function_name}()
    RETURNS TRIGGER
    LANGUAGE PLPGSQL
AS $$
{self._build_variable_declarations()}
BEGIN
{function_body}
END;
$$
"""  # noqa: S608


class PostgresWriteLogRecordingTriggerInstaller(WriteLogRecordingTriggerInstaller):
    def __init__(
        self,
        connection,
        db_metadata: DatabaseMetadataWrapper,
    ):
        super().__init__(connection=connection, db_metadata=db_metadata)

    def _install_base62_function(self, cursor: CursorWrapper):
        cursor.execute("""
CREATE OR REPLACE FUNCTION base62(n_char integer) RETURNS text AS $$
DECLARE
    alphabet text := '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
    result text := '';
    i integer;
BEGIN
    FOR i IN 1..n_char LOOP
        result := result || substr(alphabet, 1 + floor(random() * 62)::integer, 1);
    END LOOP;
    RETURN result;
END;
$$ LANGUAGE plpgsql;
""")

    @override
    def install_triggers(self, table: str, cursor: CursorWrapper):
        self._install_base62_function(cursor)

        trigger_builder = PostgresTriggerBuilder(
            table=table, db_metadata=self.db_metadata, cursor=cursor
        )

        create_trigger_function_command = trigger_builder.build()

        logger.debug(create_trigger_function_command)
        cursor.execute(create_trigger_function_command)

        for write_log_event_type in WriteLogEventTypes:
            trigger_name = (
                f"lamindb_writelog_{table}_{write_log_event_type.name.lower()}_tr"
            )

            cursor.execute(f"""
            CREATE OR REPLACE TRIGGER {trigger_name}
            AFTER {write_log_event_type.name}
            ON {table}
            FOR EACH ROW EXECUTE PROCEDURE {trigger_builder.function_name}();
            """)


def create_writelog_recording_trigger_installer(
    connection: BaseDatabaseWrapper,
) -> WriteLogRecordingTriggerInstaller:
    if connection.vendor == "postgresql":
        return PostgresWriteLogRecordingTriggerInstaller(
            connection=connection, db_metadata=PostgresDatabaseMetadataWrapper()
        )
    # TODO: add a write-log-recording trigger installer for SQLite
    else:
        raise ValueError(
            f"Write logging is not supported for vendor '{connection.vendor}'"
        )
