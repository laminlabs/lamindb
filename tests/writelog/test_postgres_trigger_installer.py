from typing import TYPE_CHECKING, Generator, cast
from unittest.mock import ANY, MagicMock

import pytest
from django.db import connection as django_connection_proxy
from django.db import transaction
from django.db.backends.utils import CursorWrapper
from lamindb.core.writelog._db_metadata_wrapper import (
    PostgresDatabaseMetadataWrapper,
)
from lamindb.core.writelog._trigger_installer import (
    FOREIGN_KEYS_LIST_COLUMN_NAME,
    PostgresWriteLogDBRecordingTriggerInstaller,
    WriteLogEventTypes,
)
from lamindb.core.writelog._types import TableUID, UIDColumns
from lamindb.models.record import Space
from lamindb.models.writelog import WriteLog, WriteLogMigrationState, WriteLogTableState
from typing_extensions import override

if TYPE_CHECKING:
    from django.db.backends.base.base import BaseDatabaseWrapper

django_connection = cast("BaseDatabaseWrapper", django_connection_proxy)


class FakeMetadataWrapper(PostgresDatabaseMetadataWrapper):
    """A fake DB metadata wrapper that allows us to control which database tables the installer will see and target."""

    def __init__(self):
        super().__init__()
        self._tables_with_triggers = set()
        self._db_tables = set()
        self._many_to_many_tables = set()
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


def fetch_row_id_by_uid(
    table: str,
    id_columns: list[str],
    uid_columns: dict[str, str],
    cursor: CursorWrapper,
) -> list[int]:
    where_clause = " AND ".join(f"{k} = '{v}'" for (k, v) in uid_columns.items())

    cursor.execute(f"SELECT {', '.join(id_columns)} FROM {table} WHERE {where_clause}")  # noqa: S608
    row = cursor.fetchone()
    assert row is not None

    return [int(r) for r in row]


@pytest.fixture(scope="function", autouse=True)
def write_log_state():
    # Clean up after any tests that modify write logs
    WriteLog.objects.all().delete()
    WriteLogTableState.objects.all().delete()
    WriteLogMigrationState.objects.all().delete()

    yield

    WriteLog.objects.all().delete()
    WriteLogTableState.objects.all().delete()
    WriteLogMigrationState.objects.all().delete()


@pytest.fixture(scope="function")
def table_a(write_log_state) -> Generator[str, None, None]:
    table_name = "write_log_test_table_a"
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS {table_name}
(id SERIAL PRIMARY KEY, str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield table_name

    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")


@pytest.fixture(scope="function")
def table_b(write_log_state) -> Generator[str, None, None]:
    table_name = "write_log_test_table_b"
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS {table_name}
(id SERIAL PRIMARY KEY, str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield table_name

    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")


@pytest.fixture(scope="function")
def table_c(write_log_state) -> Generator[str, None, None]:
    table_name = "write_log_test_table_c"
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS {table_name}
(id SERIAL PRIMARY KEY, str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield table_name

    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")


@pytest.mark.pg_integration
def test_updating_write_log_triggers_installs_table_state(table_a, table_b, table_c):
    assert len(set(WriteLogTableState.objects.all())) == 0

    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b, table_c})
    fake_db_metadata.set_tables_with_installed_triggers({table_a, table_c})

    installer = PostgresWriteLogDBRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )
    installer.install_triggers = MagicMock(return_value=None)

    installer.update_write_log_triggers()

    new_table_states = set(WriteLogTableState.objects.all())

    assert len(new_table_states) == 3
    assert {ts.table_name for ts in new_table_states} == {
        table_a,
        table_b,
        table_c,
    }
    assert all(ts.backfilled is False for ts in new_table_states)


@pytest.mark.pg_integration
def test_update_write_log_triggers_only_install_table_state_once(
    table_a, table_b, table_c
):
    assert len(set(WriteLogTableState.objects.all())) == 0

    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b, table_c})
    fake_db_metadata.set_tables_with_installed_triggers({table_a, table_c})

    installer = PostgresWriteLogDBRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )
    installer.install_triggers = MagicMock(return_value=None)

    # Run update_write_log_triggers() twice.
    installer.update_write_log_triggers(update_all=True)
    installer.update_write_log_triggers(update_all=True)

    new_table_states = set(WriteLogTableState.objects.all())

    # There should only be one table state per table, even though we've run more than once.
    assert len(new_table_states) == 3
    assert {ts.table_name for ts in new_table_states} == {
        table_a,
        table_b,
        table_c,
    }
    assert all(ts.backfilled is False for ts in new_table_states)


@pytest.mark.pg_integration
def test_update_triggers_installs_migration_state(table_a, table_b, table_c):
    assert len(set(WriteLogMigrationState.objects.all())) == 0

    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b, table_c})
    fake_db_metadata.set_tables_with_installed_triggers(set())

    installer = PostgresWriteLogDBRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )

    installer.install_triggers = MagicMock(return_value=None)

    installer.update_write_log_triggers()

    new_migration_states = set(WriteLogMigrationState.objects.all())

    assert len(new_migration_states) == 1

    # Updating triggers a second time shouldn't add more migration state, since
    # no migrations have been performed between the two updates.

    installer.update_write_log_triggers()

    new_migration_states = set(WriteLogMigrationState.objects.all())

    assert len(new_migration_states) == 1


@pytest.mark.pg_integration
def test_update_write_log_triggers_skips_existing_triggers(table_a, table_b, table_c):
    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b, table_c})
    fake_db_metadata.set_tables_with_installed_triggers({table_a, table_c})

    installer = PostgresWriteLogDBRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )

    installer.install_triggers = MagicMock(return_value=None)

    installer.update_write_log_triggers()

    installer.install_triggers.assert_called_once_with(table_b, ANY)


@pytest.fixture(scope="function")
def no_uid_pg_table(write_log_state) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS write_log_table_test_no_uid
(id SERIAL PRIMARY KEY, str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield "write_log_table_test_no_uid"

    cursor.execute("DROP TABLE IF EXISTS write_log_table_test_no_uid")


@pytest.fixture(scope="function")
def foreign_key_to_no_uid_table(
    write_log_state, no_uid_pg_table
) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS write_log_table_test_no_uid_fk
(id SERIAL PRIMARY KEY, primary_id INT, CONSTRAINT fk FOREIGN KEY (primary_id) REFERENCES {no_uid_pg_table}(id));
""")

    yield "write_log_table_test_no_uid_fk"

    cursor.execute("DROP TABLE IF EXISTS write_log_table_test_no_uid_fk")


@pytest.mark.pg_integration
def test_foreign_key_to_table_without_uid_fails(
    no_uid_pg_table, foreign_key_to_no_uid_table
):
    with pytest.raises(ValueError):
        _update_write_log_triggers({no_uid_pg_table, foreign_key_to_no_uid_table})


@pytest.mark.pg_integration
def test_sql_injectable_table_names_fail(table_a, table_b):
    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b})
    fake_db_metadata.set_tables_with_installed_triggers(set())

    installer = PostgresWriteLogDBRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )

    with pytest.raises(ValueError):
        installer.install_triggers(
            "bad_tabl; DROP TABLE foo", django_connection.cursor()
        )


@pytest.fixture(scope="function")
def simple_pg_table(write_log_state) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS write_log_table_test_a
(id SERIAL PRIMARY KEY, uid VARCHAR(8), str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield "write_log_table_test_a"

    cursor.execute("DROP TABLE IF EXISTS write_log_table_test_a")


def _update_write_log_triggers(
    table_list: set[str],
    many_to_many_tables: set[str] | None = None,
    uid_columns: dict[str, UIDColumns] | None = None,
):
    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables(table_list)
    fake_db_metadata.set_tables_with_installed_triggers(set())

    if many_to_many_tables is not None:
        fake_db_metadata.set_many_to_many_db_tables(many_to_many_tables)

    if uid_columns is not None:
        for table, uid_columns_for_table in uid_columns.items():
            fake_db_metadata.set_uid_columns(table, uid_columns_for_table)

    installer = PostgresWriteLogDBRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )

    installer.update_write_log_triggers()


@pytest.mark.pg_integration
def test_simple_table_trigger(simple_pg_table: str):
    _update_write_log_triggers({simple_pg_table})

    cursor = django_connection.cursor()

    cursor.execute("SELECT * FROM pg_trigger WHERE tgname LIKE 'lamindb_writelog_%';")

    triggers = cursor.fetchall()

    assert len(triggers) == 3

    cursor.execute(
        f"INSERT INTO {simple_pg_table} (uid, str_col, bool_col, int_col) "  # noqa: S608
        "VALUES ('abc123', 'hello world', false, 42)"
    )
    cursor.execute(
        f"INSERT INTO {simple_pg_table} (uid, str_col, bool_col, int_col) "  # noqa: S608
        "VALUES ('def456', 'another row', true, 99)"
    )
    cursor.execute(f"UPDATE {simple_pg_table} SET int_col=22 WHERE uid = 'def456'")  # noqa: S608
    cursor.execute(f"DELETE FROM {simple_pg_table} WHERE uid = 'abc123'")  # noqa: S608

    write_log = WriteLog.objects.all().order_by("seqno")

    assert len(write_log) == 4

    assert [h.table.table_name == simple_pg_table for h in write_log]

    assert write_log[0].record_uid == ["abc123"]
    assert write_log[0].record_data == {
        "int_col": 42,
        "str_col": "hello world",
        "bool_col": False,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[0].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[1].record_uid == ["def456"]
    assert write_log[1].record_data == {
        "int_col": 99,
        "str_col": "another row",
        "bool_col": True,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[1].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[2].record_uid == ["def456"]
    assert write_log[2].record_data == {
        "int_col": 22,
        "str_col": "another row",
        "bool_col": True,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[2].event_type == WriteLogEventTypes.UPDATE.value

    assert write_log[3].record_uid == ["abc123"]
    assert write_log[3].record_data is None
    assert write_log[3].event_type == WriteLogEventTypes.DELETE.value


@pytest.fixture(scope="function")
def foreignkey_pg_table(simple_pg_table) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS write_log_table_test_b
(
    id SERIAL PRIMARY KEY, uid VARCHAR(8), table_a_id int,
    CONSTRAINT fk_table_a FOREIGN KEY(table_a_id) REFERENCES {simple_pg_table}(id)
    ON DELETE CASCADE
);
""")

    yield "write_log_table_test_b"

    cursor.execute("DROP TABLE IF EXISTS write_log_table_test_b")


@pytest.mark.pg_integration
def test_triggers_with_foreign_keys(simple_pg_table, foreignkey_pg_table):
    _update_write_log_triggers({simple_pg_table, foreignkey_pg_table})

    cursor = django_connection.cursor()

    with transaction.atomic():
        cursor.execute(
            f"INSERT INTO {simple_pg_table} (uid, str_col, bool_col, int_col) "  # noqa: S608
            "VALUES ('abc123', 'hello world', false, 42)"
        )

        table_a_row_id = fetch_row_id_by_uid(
            table=simple_pg_table,
            id_columns=["id"],
            uid_columns={"uid": "abc123"},
            cursor=cursor,
        )[0]

        cursor.execute(
            f"INSERT INTO {foreignkey_pg_table} (uid, table_a_id) VALUES ('foo333', {table_a_row_id})"  # noqa: S608
        )
    cursor.execute(
        f"INSERT INTO {foreignkey_pg_table} (uid, table_a_id) VALUES ('bar444', NULL)"  # noqa: S608
    )
    cursor.execute(f"DELETE FROM {simple_pg_table} WHERE uid = 'abc123'")  # noqa: S608

    write_log = WriteLog.objects.all().order_by("seqno")

    assert len(write_log) == 5

    assert write_log[0].table.table_name == simple_pg_table
    assert write_log[0].record_uid == ["abc123"]
    assert write_log[0].record_data == {
        "int_col": 42,
        "str_col": "hello world",
        "bool_col": False,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[0].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[1].table.table_name == foreignkey_pg_table
    assert write_log[1].record_uid == ["foo333"]
    assert write_log[1].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["table_a_id"], {"uid": "abc123"}]
        ],
    }
    assert write_log[1].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[2].table.table_name == foreignkey_pg_table
    assert write_log[2].record_uid == ["bar444"]
    assert write_log[2].event_type == WriteLogEventTypes.INSERT.value
    assert write_log[2].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["table_a_id"], {"uid": None}]
        ],
    }

    assert write_log[3].table.table_name == simple_pg_table
    assert write_log[3].record_uid == ["abc123"]
    assert write_log[3].event_type == WriteLogEventTypes.DELETE.value
    assert write_log[3].record_data is None

    assert write_log[4].table.table_name == foreignkey_pg_table
    assert write_log[4].record_uid == ["foo333"]
    assert write_log[4].event_type == WriteLogEventTypes.DELETE.value
    assert write_log[4].record_data is None


@pytest.fixture(scope="function")
def self_referential_pg_table() -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS write_log_table_test_c
(
    id SERIAL PRIMARY KEY, uid VARCHAR(8), parent_id int,
    CONSTRAINT fk_parent FOREIGN KEY(parent_id) REFERENCES write_log_table_test_c(id)
    ON DELETE CASCADE
);
""")

    yield "write_log_table_test_c"

    cursor.execute("DROP TABLE IF EXISTS write_log_table_test_c")


@pytest.mark.pg_integration
def test_triggers_with_self_references(self_referential_pg_table):
    _update_write_log_triggers({self_referential_pg_table})

    cursor = django_connection.cursor()

    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) "  # noqa: S608
        "VALUES ('abc123', NULL)"
    )
    parent_row_id = fetch_row_id_by_uid(
        table=self_referential_pg_table,
        id_columns=["id"],
        uid_columns={"uid": "abc123"},
        cursor=cursor,
    )[0]

    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) VALUES ('def345', {parent_row_id})"  # noqa: S608
    )

    write_log = WriteLog.objects.all().order_by("seqno")

    assert len(write_log) == 2

    assert write_log[0].table.table_name == self_referential_pg_table
    assert write_log[0].record_uid == ["abc123"]
    assert write_log[0].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["parent_id"], {"uid": None}]
        ],
    }
    assert write_log[0].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[1].table.table_name == self_referential_pg_table
    assert write_log[1].record_uid == ["def345"]
    assert write_log[1].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["parent_id"], {"uid": "abc123"}]
        ],
    }
    assert write_log[1].event_type == WriteLogEventTypes.INSERT.value


@pytest.fixture(scope="function")
def composite_primary_key_table() -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS write_log_table_test_d
(
    key_a INT,
    key_b INT,
    uid VARCHAR(8),
    PRIMARY KEY (key_a, key_b)
);
""")

    yield "write_log_table_test_d"

    cursor.execute("DROP TABLE IF EXISTS write_log_table_test_d")


@pytest.fixture(scope="function")
def composite_foreign_key_table(
    composite_primary_key_table,
) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS write_log_table_test_e
(
    id SERIAL PRIMARY KEY,
    table_d_key_a INT,
    table_d_key_b INT,
    uid VARCHAR(8),
    CONSTRAINT fk_parent FOREIGN KEY(table_d_key_a, table_d_key_b) REFERENCES {composite_primary_key_table} (key_a, key_b)
);
""")

    yield "write_log_table_test_e"

    cursor.execute("DROP TABLE IF EXISTS write_log_table_test_e")


@pytest.mark.pg_integration
def test_triggers_with_composite_primary_key(
    composite_primary_key_table, composite_foreign_key_table
):
    _update_write_log_triggers(
        {composite_primary_key_table, composite_foreign_key_table}
    )

    cursor = django_connection.cursor()

    cursor.execute(
        f"INSERT INTO {composite_primary_key_table} (key_a, key_b, uid) "  # noqa: S608
        "VALUES (42, 21, 'xyz123')"
    )
    cursor.execute(
        f"INSERT INTO {composite_foreign_key_table} (table_d_key_a, table_d_key_b, uid) "  # noqa: S608
        f"VALUES (42, 21, 'qrs234')"
    )

    write_log = WriteLog.objects.all().order_by("seqno")

    assert len(write_log) == 2

    assert write_log[0].table.table_name == composite_primary_key_table
    assert write_log[0].record_uid == ["xyz123"]
    assert write_log[0].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[0].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[1].table.table_name == composite_foreign_key_table
    assert write_log[1].record_uid == ["qrs234"]
    assert write_log[1].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [
                write_log[0].table.id,
                ["table_d_key_a", "table_d_key_b"],
                {"uid": "xyz123"},
            ]
        ],
    }
    assert write_log[1].event_type == WriteLogEventTypes.INSERT.value


@pytest.fixture(scope="function")
def many_to_many_table(
    composite_primary_key_table, self_referential_pg_table
) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS write_log_table_many_to_many_test
(
    id SERIAL PRIMARY KEY,
    composite_key_a INT,
    composite_key_b INT,
    self_ref_id INT,
    CONSTRAINT fk_one FOREIGN KEY(composite_key_a, composite_key_b) REFERENCES {composite_primary_key_table} (key_a, key_b),
    CONSTRAINT fk_two FOREIGN KEY(self_ref_id) REFERENCES {self_referential_pg_table} (id)
);
""")

    yield "write_log_table_many_to_many_test"

    cursor.execute("DROP TABLE IF EXISTS write_log_table_many_to_many_test")


@pytest.mark.pg_integration
def test_triggers_with_many_to_many_tables(
    many_to_many_table, composite_primary_key_table, self_referential_pg_table
):
    _update_write_log_triggers(
        table_list={
            composite_primary_key_table,
            self_referential_pg_table,
            many_to_many_table,
        },
        many_to_many_tables={many_to_many_table},
    )

    cursor = django_connection.cursor()

    cursor.execute(
        f"INSERT INTO {composite_primary_key_table} (key_a, key_b, uid) "  # noqa: S608
        "VALUES (42, 21, 'xyz123')"
    )

    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) "  # noqa: S608
        "VALUES ('abc123', NULL)"
    )

    composite_row_id = fetch_row_id_by_uid(
        table=composite_primary_key_table,
        id_columns=["key_a", "key_b"],
        uid_columns={"uid": "xyz123"},
        cursor=cursor,
    )
    self_ref_row_id = fetch_row_id_by_uid(
        table=self_referential_pg_table,
        id_columns=["id"],
        uid_columns={"uid": "abc123"},
        cursor=cursor,
    )[0]

    cursor.execute(
        f"INSERT INTO {many_to_many_table} (composite_key_a, composite_key_b, self_ref_id) VALUES ({composite_row_id[0]}, {composite_row_id[1]}, {self_ref_row_id})"  # noqa: S608
    )

    cursor.execute(
        f"DELETE FROM {many_to_many_table} WHERE self_ref_id={self_ref_row_id}"  # noqa: S608
    )

    cursor.execute(
        f"DELETE FROM {self_referential_pg_table} WHERE id={self_ref_row_id}"  # noqa: S608
    )

    write_log = WriteLog.objects.all().order_by("seqno")

    assert len(write_log) == 5

    assert write_log[0].table.table_name == composite_primary_key_table
    assert write_log[0].record_uid == ["xyz123"]
    assert write_log[0].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[0].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[1].table.table_name == self_referential_pg_table
    assert write_log[1].record_uid == ["abc123"]
    assert write_log[1].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[1].table.id, ["parent_id"], {"uid": None}]
        ]
    }
    assert write_log[1].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[2].table.table_name == many_to_many_table
    assert write_log[2].record_uid == [
        [
            write_log[0].table.id,
            ["composite_key_a", "composite_key_b"],
            {"uid": "xyz123"},
        ],
        [write_log[1].table.id, ["self_ref_id"], {"uid": "abc123"}],
    ]
    assert write_log[2].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [
                write_log[0].table.id,
                ["composite_key_a", "composite_key_b"],
                {"uid": "xyz123"},
            ],
            [write_log[1].table.id, ["self_ref_id"], {"uid": "abc123"}],
        ]
    }
    assert write_log[2].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[3].table.table_name == many_to_many_table
    assert write_log[3].record_uid == [
        [
            write_log[0].table.id,
            ["composite_key_a", "composite_key_b"],
            {"uid": "xyz123"},
        ],
        [write_log[1].table.id, ["self_ref_id"], {"uid": "abc123"}],
    ]
    assert write_log[3].record_data is None
    assert write_log[3].event_type == WriteLogEventTypes.DELETE.value

    assert write_log[4].table.table_name == self_referential_pg_table
    assert write_log[4].record_uid == ["abc123"]
    assert write_log[4].record_data is None
    assert write_log[4].event_type == WriteLogEventTypes.DELETE.value


@pytest.fixture(scope="function")
def compound_uid_table() -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS compound_uid_test
(
    id SERIAL PRIMARY KEY,
    uid_1 VARCHAR(8),
    uid_2 VARCHAR(8)
);
""")

    yield "compound_uid_test"

    cursor.execute("DROP TABLE IF EXISTS compound_uid_test")


@pytest.fixture(scope="function")
def compound_uid_child(compound_uid_table) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS compound_uid_child_test
(
    id SERIAL PRIMARY KEY,
    uid VARCHAR(8),
    parent_id INT,
    CONSTRAINT fk FOREIGN KEY(parent_id) REFERENCES {compound_uid_table} (id)
);
""")

    yield "compound_uid_child_test"

    cursor.execute("DROP TABLE IF EXISTS compound_uid_child_test")


@pytest.mark.pg_integration
def test_triggers_with_compound_table_uid(compound_uid_table, compound_uid_child):
    _update_write_log_triggers(
        table_list={compound_uid_table, compound_uid_child},
        uid_columns={
            compound_uid_table: [
                TableUID(
                    source_table_name=compound_uid_table,
                    uid_columns=["uid_1", "uid_2"],
                    key_constraint=None,
                )
            ]
        },
    )

    cursor = django_connection.cursor()

    cursor.execute(
        f"INSERT INTO {compound_uid_table} (uid_1, uid_2) "  # noqa: S608
        "VALUES ('badf00d', 'dadb0d')"
    )

    row_id = fetch_row_id_by_uid(
        table=compound_uid_table,
        id_columns=["id"],
        uid_columns={"uid_1": "badf00d", "uid_2": "dadb0d"},
        cursor=cursor,
    )[0]

    cursor.execute(
        f"INSERT INTO {compound_uid_child} (uid, parent_id) VALUES ('abc123', {row_id})"  # noqa: S608
    )

    cursor.execute(
        f"INSERT INTO {compound_uid_child} (uid, parent_id) VALUES ('def456', NULL)"  # noqa: S608
    )

    cursor.execute(f"DELETE FROM {compound_uid_child} WHERE uid = 'abc123'")  # noqa: S608
    cursor.execute(
        f"DELETE FROM {compound_uid_table} WHERE uid_1 = 'badf00d' AND uid_2 = 'dadb0d'"  # noqa: S608
    )

    write_log = WriteLog.objects.all().order_by("seqno")

    assert len(write_log) == 5

    assert write_log[0].table.table_name == compound_uid_table
    assert write_log[0].record_uid == ["badf00d", "dadb0d"]
    assert write_log[0].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[0].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[1].table.table_name == compound_uid_child
    assert write_log[1].record_uid == ["abc123"]
    assert write_log[1].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [
                write_log[0].table.id,
                ["parent_id"],
                {"uid_1": "badf00d", "uid_2": "dadb0d"},
            ]
        ],
    }
    assert write_log[1].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[2].table.table_name == compound_uid_child
    assert write_log[2].record_uid == ["def456"]
    assert write_log[2].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["parent_id"], {"uid_1": None, "uid_2": None}]
        ],
    }
    assert write_log[2].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[3].table.table_name == compound_uid_child
    assert write_log[3].record_uid == ["abc123"]
    assert write_log[3].record_data is None
    assert write_log[3].event_type == WriteLogEventTypes.DELETE.value

    assert write_log[4].table.table_name == compound_uid_table
    assert write_log[4].record_uid == ["badf00d", "dadb0d"]
    assert write_log[4].record_data is None
    assert write_log[4].event_type == WriteLogEventTypes.DELETE.value


@pytest.fixture(scope="function")
def compound_uid_many_to_many(compound_uid_table) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS compound_uid_many_to_many_test
(
    id SERIAL PRIMARY KEY,
    table_a_id INT,
    table_b_id INT,
    CONSTRAINT fk_a FOREIGN KEY(table_a_id) REFERENCES {compound_uid_table} (id),
    CONSTRAINT fk_b FOREIGN KEY(table_b_id) REFERENCES {compound_uid_table} (id)

);
""")

    yield "compound_uid_many_to_many_test"

    cursor.execute("DROP TABLE IF EXISTS compound_uid_many_to_many_test")


@pytest.mark.pg_integration
def test_triggers_many_to_many_to_compound_uid_with_self_links(
    compound_uid_table, compound_uid_many_to_many
):
    _update_write_log_triggers(
        table_list={compound_uid_table, compound_uid_many_to_many},
        many_to_many_tables={compound_uid_many_to_many},
        uid_columns={
            compound_uid_table: [
                TableUID(
                    source_table_name=compound_uid_table,
                    uid_columns=["uid_1", "uid_2"],
                    key_constraint=None,
                )
            ]
        },
    )

    cursor = django_connection.cursor()

    cursor.execute(
        f"INSERT INTO {compound_uid_table} (uid_1, uid_2) "  # noqa: S608
        "VALUES ('rec1_1', 'rec1_2')"
    )

    cursor.execute(
        f"INSERT INTO {compound_uid_table} (uid_1, uid_2) "  # noqa: S608
        "VALUES ('rec2_1', 'rec2_2')"
    )

    rec_1_id = fetch_row_id_by_uid(
        table=compound_uid_table,
        id_columns=["id"],
        uid_columns={"uid_1": "rec1_1", "uid_2": "rec1_2"},
        cursor=cursor,
    )[0]
    rec_2_id = fetch_row_id_by_uid(
        table=compound_uid_table,
        id_columns=["id"],
        uid_columns={"uid_1": "rec2_1", "uid_2": "rec2_2"},
        cursor=cursor,
    )[0]

    cursor.execute(
        f"INSERT INTO {compound_uid_many_to_many} (table_a_id, table_b_id) VALUES ({rec_1_id}, {rec_2_id})"  # noqa: S608
    )
    cursor.execute(
        f"INSERT INTO {compound_uid_many_to_many} (table_a_id, table_b_id) VALUES (NULL, {rec_1_id})"  # noqa: S608
    )
    cursor.execute(
        f"INSERT INTO {compound_uid_many_to_many} (table_a_id, table_b_id) VALUES ({rec_2_id}, NULL)"  # noqa: S608
    )
    cursor.execute(
        f"INSERT INTO {compound_uid_many_to_many} (table_a_id, table_b_id) VALUES (NULL, NULL)"  # noqa: S608
    )

    cursor.execute(
        f"DELETE FROM {compound_uid_many_to_many} WHERE table_a_id = {rec_1_id} AND table_b_id = {rec_2_id}"  # noqa: S608
    )

    write_log = WriteLog.objects.all().order_by("seqno")

    assert len(write_log) == 7

    assert write_log[0].table.table_name == compound_uid_table
    assert write_log[0].record_uid == ["rec1_1", "rec1_2"]
    assert write_log[0].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[0].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[1].table.table_name == compound_uid_table
    assert write_log[1].record_uid == ["rec2_1", "rec2_2"]
    assert write_log[1].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[1].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[2].table.table_name == compound_uid_many_to_many
    assert write_log[2].record_uid == [
        [write_log[0].table.id, ["table_a_id"], {"uid_1": "rec1_1", "uid_2": "rec1_2"}],
        [write_log[0].table.id, ["table_b_id"], {"uid_1": "rec2_1", "uid_2": "rec2_2"}],
    ]
    assert write_log[2].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [
                write_log[0].table.id,
                ["table_a_id"],
                {"uid_1": "rec1_1", "uid_2": "rec1_2"},
            ],
            [
                write_log[0].table.id,
                ["table_b_id"],
                {"uid_1": "rec2_1", "uid_2": "rec2_2"},
            ],
        ]
    }
    assert write_log[2].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[3].table.table_name == compound_uid_many_to_many
    assert write_log[3].record_uid == [
        [write_log[0].table.id, ["table_a_id"], {"uid_1": None, "uid_2": None}],
        [write_log[0].table.id, ["table_b_id"], {"uid_1": "rec1_1", "uid_2": "rec1_2"}],
    ]
    assert write_log[3].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["table_a_id"], {"uid_1": None, "uid_2": None}],
            [
                write_log[0].table.id,
                ["table_b_id"],
                {"uid_1": "rec1_1", "uid_2": "rec1_2"},
            ],
        ]
    }
    assert write_log[3].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[4].table.table_name == compound_uid_many_to_many
    assert write_log[4].record_uid == [
        [write_log[0].table.id, ["table_a_id"], {"uid_1": "rec2_1", "uid_2": "rec2_2"}],
        [write_log[0].table.id, ["table_b_id"], {"uid_1": None, "uid_2": None}],
    ]
    assert write_log[4].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [
                write_log[0].table.id,
                ["table_a_id"],
                {"uid_1": "rec2_1", "uid_2": "rec2_2"},
            ],
            [write_log[0].table.id, ["table_b_id"], {"uid_1": None, "uid_2": None}],
        ]
    }
    assert write_log[4].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[5].table.table_name == compound_uid_many_to_many
    assert write_log[5].record_uid == [
        [write_log[0].table.id, ["table_a_id"], {"uid_1": None, "uid_2": None}],
        [write_log[0].table.id, ["table_b_id"], {"uid_1": None, "uid_2": None}],
    ]
    assert write_log[5].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["table_a_id"], {"uid_1": None, "uid_2": None}],
            [write_log[0].table.id, ["table_b_id"], {"uid_1": None, "uid_2": None}],
        ]
    }
    assert write_log[5].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[6].table.table_name == compound_uid_many_to_many
    assert write_log[6].record_uid == [
        [write_log[0].table.id, ["table_a_id"], {"uid_1": "rec1_1", "uid_2": "rec1_2"}],
        [write_log[0].table.id, ["table_b_id"], {"uid_1": "rec2_1", "uid_2": "rec2_2"}],
    ]
    assert write_log[6].record_data is None
    assert write_log[6].event_type == WriteLogEventTypes.DELETE.value


@pytest.fixture(scope="function")
def fake_space():
    space = Space(name="my_fake_space", uid="fakespace").save()  # type: ignore

    yield space

    space.delete()


@pytest.fixture(scope="function")
def table_with_space_ref(fake_space):
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS write_log_space_ref_test
(
    id SERIAL PRIMARY KEY,
    uid VARCHAR(8),
    space_id INT,
    CONSTRAINT fk FOREIGN KEY(space_id) REFERENCES lamindb_space (id)
);
""")

    yield "write_log_space_ref_test"

    cursor.execute("DROP TABLE IF EXISTS write_log_space_ref_test")


@pytest.mark.pg_integration
def test_write_log_records_space_uids_properly(table_with_space_ref, fake_space):
    cursor = django_connection.cursor()

    _update_write_log_triggers(
        table_list={table_with_space_ref},
    )

    cursor = django_connection.cursor()

    cursor.execute(
        f"INSERT INTO {table_with_space_ref} (uid, space_id) "  # noqa: S608
        f"VALUES ('A', {fake_space.id})"
    )

    cursor.execute(
        f"INSERT INTO {table_with_space_ref} (uid, space_id) "  # noqa: S608
        f"VALUES ('B', NULL)"
    )

    write_log = WriteLog.objects.all().order_by("seqno")

    assert len(write_log) == 2

    assert write_log[0].table.table_name == table_with_space_ref
    assert write_log[0].record_uid == ["A"]
    assert write_log[0].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[0].event_type == WriteLogEventTypes.INSERT.value
    assert write_log[0].space_uid == "fakespace"

    assert write_log[1].table.table_name == table_with_space_ref
    assert write_log[1].record_uid == ["B"]
    assert write_log[1].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert write_log[1].event_type == WriteLogEventTypes.INSERT.value
    assert write_log[1].space_uid is None


@pytest.mark.pg_integration
def test_write_log_install_triggers_on_existing_lamindb_models():
    cursor = django_connection.cursor()

    try:
        installer = PostgresWriteLogDBRecordingTriggerInstaller(
            connection=django_connection, db_metadata=PostgresDatabaseMetadataWrapper()
        )
        installer.update_write_log_triggers()
    finally:
        # Drop all write log triggers
        cursor.execute("""
DO $$
DECLARE
    r RECORD;
BEGIN
    FOR r IN
        SELECT
            tgname AS trigger_name,
            relname AS table_name,
            nspname AS schema_name
        FROM pg_trigger t
        JOIN pg_class c ON t.tgrelid = c.oid
        JOIN pg_namespace n ON c.relnamespace = n.oid
        WHERE tgname LIKE 'lamindb_writelog_%'
          AND NOT tgisinternal
    LOOP
        EXECUTE format('DROP TRIGGER %I ON %I.%I;',
                      r.trigger_name, r.schema_name, r.table_name);
        RAISE NOTICE 'Dropped trigger % on %.%',
                    r.trigger_name, r.schema_name, r.table_name;
    END LOOP;
END;
$$;
                       """)
