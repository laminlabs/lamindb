import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Generator, cast
from unittest.mock import ANY, MagicMock

import pytest
from django.db import connection as django_connection_proxy
from django.db import transaction
from django.db.backends.utils import CursorWrapper
from lamindb.core.writelog._constants import FOREIGN_KEYS_LIST_COLUMN_NAME
from lamindb.core.writelog._db_metadata_wrapper import (
    PostgresDatabaseMetadataWrapper,
)
from lamindb.core.writelog._trigger_installer import (
    PostgresWriteLogRecordingTriggerInstaller,
    WriteLogEventTypes,
)
from lamindb.core.writelog._types import Column, ColumnType, TableUID, UIDColumns
from lamindb.models.artifact import Artifact
from lamindb.models.run import Run
from lamindb.models.sqlrecord import Space
from lamindb.models.transform import Transform
from lamindb.models.writelog import (
    WriteLog,
    WriteLogLock,
    WriteLogMigrationState,
    WriteLogTableState,
)
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

    installer = PostgresWriteLogRecordingTriggerInstaller(
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

    assert not WriteLogTableState.objects.get(table_name=table_a).backfilled
    assert WriteLogTableState.objects.get(table_name=table_b).backfilled
    assert not WriteLogTableState.objects.get(table_name=table_c).backfilled


@pytest.mark.pg_integration
def test_update_write_log_triggers_only_install_table_state_once(
    table_a, table_b, table_c
):
    assert len(set(WriteLogTableState.objects.all())) == 0

    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b, table_c})
    fake_db_metadata.set_tables_with_installed_triggers({table_a, table_c})

    installer = PostgresWriteLogRecordingTriggerInstaller(
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

    assert all(t.backfilled for t in new_table_states)


@pytest.mark.pg_integration
def test_update_triggers_installs_migration_state(table_a, table_b, table_c):
    assert len(set(WriteLogMigrationState.objects.all())) == 0

    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b, table_c})
    fake_db_metadata.set_tables_with_installed_triggers(set())

    installer = PostgresWriteLogRecordingTriggerInstaller(
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

    installer = PostgresWriteLogRecordingTriggerInstaller(
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

    installer = PostgresWriteLogRecordingTriggerInstaller(
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
    tables_with_installed_triggers: set[str] | None = None,
):
    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables(table_list)

    if tables_with_installed_triggers is not None:
        fake_db_metadata.set_tables_with_installed_triggers(
            tables_with_installed_triggers
        )

    if many_to_many_tables is not None:
        fake_db_metadata.set_many_to_many_db_tables(many_to_many_tables)

    if uid_columns is not None:
        for table, uid_columns_for_table in uid_columns.items():
            fake_db_metadata.set_uid_columns(table, uid_columns_for_table)

    installer = PostgresWriteLogRecordingTriggerInstaller(
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

    write_log = WriteLog.objects.all().order_by("id")

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
def write_log_lock_locked():
    write_log_lock = WriteLogLock.load()

    assert write_log_lock is not None

    write_log_lock.lock()

    yield write_log_lock

    write_log_lock.unlock()


@pytest.mark.pg_integration
def test_trigger_doesnt_fire_when_write_lock_is_locked(
    simple_pg_table: str, write_log_lock_locked
):
    _update_write_log_triggers({simple_pg_table})

    cursor = django_connection.cursor()

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

    write_log = WriteLog.objects.all().order_by("id")

    assert len(write_log) == 0


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

    write_log = WriteLog.objects.all().order_by("id")

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

    write_log = WriteLog.objects.all().order_by("id")

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

    write_log = WriteLog.objects.all().order_by("id")

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

    write_log = WriteLog.objects.all().order_by("id")

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
                    uid_columns=[
                        Column(name="uid_1", type=ColumnType.STR, ordinal_position=1),
                        Column(name="uid_2", type=ColumnType.STR, ordinal_position=3),
                    ],
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

    write_log = WriteLog.objects.all().order_by("id")

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
                    uid_columns=[
                        Column(name="uid_1", type=ColumnType.STR, ordinal_position=1),
                        Column(name="uid_2", type=ColumnType.STR, ordinal_position=3),
                    ],
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

    write_log = WriteLog.objects.all().order_by("id")

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


@pytest.fixture(scope="function")
def backfill_table_a(write_log_state):
    table_name = "backfill_table_a"

    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS {table_name}
(
    id SERIAL PRIMARY KEY,
    uid VARCHAR(20),
    data VARCHAR(100)
);
""")

    yield table_name

    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")


@pytest.fixture(scope="function")
def backfill_table_b(write_log_state, backfill_table_a):
    table_name = "backfill_table_b"

    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS {table_name}
(
    id SERIAL PRIMARY KEY,
    uid VARCHAR(20),
    table_a_id INT,
    data VARCHAR(100),
    CONSTRAINT fk FOREIGN KEY (table_a_id) REFERENCES {backfill_table_a}(id)
);
""")

    yield table_name

    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")


@pytest.mark.pg_integration
def test_simple_backfill(backfill_table_a, backfill_table_b):
    cursor = django_connection.cursor()

    # Populate the tables with some pre-existing data.
    cursor.execute(
        f"INSERT INTO {backfill_table_a} (uid, data) VALUES ('badf00d1234', 'mydata')"  # noqa: S608
    )
    cursor.execute(
        f"INSERT INTO {backfill_table_a} (uid, data) VALUES ('deadbeef456', 'moredata')"  # noqa: S608
    )

    table_a_id = fetch_row_id_by_uid(
        backfill_table_a, ["id"], {"uid": "badf00d1234"}, cursor
    )

    cursor.execute(
        f"INSERT INTO {backfill_table_b} (uid, table_a_id, data) VALUES ('m00m00m00', {table_a_id[0]}, 'yetmoredata')"  # noqa: S608
    )

    _update_write_log_triggers(table_list={backfill_table_a, backfill_table_b})

    write_log = WriteLog.objects.all().order_by("id")

    assert len(write_log) == 3

    # backfill_table_a should be backfilled before backfill_table_b,
    # since b has a foreign key to a

    # backfill_table_a's records can be backfilled in any order
    for i in [0, 1]:
        assert write_log[i].table.table_name == backfill_table_a
        assert write_log[i].record_uid in (["badf00d1234"], ["deadbeef456"])
        assert write_log[i].record_data == {
            "data": "mydata"
            if write_log[i].record_uid == ["badf00d1234"]
            else "moredata",
            FOREIGN_KEYS_LIST_COLUMN_NAME: [],
        }
        assert write_log[i].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[2].table.table_name == backfill_table_b
    assert write_log[2].record_uid == ["m00m00m00"]
    assert write_log[2].record_data == {
        "data": "yetmoredata",
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["table_a_id"], {"uid": "badf00d1234"}]
        ],
    }
    assert write_log[2].event_type == WriteLogEventTypes.INSERT.value


@pytest.mark.pg_integration
def test_self_referential_backfill(self_referential_pg_table):
    cursor = django_connection.cursor()

    # Populate the tables with some pre-existing data.
    # Inter-record dependencies are as follows (id -> parent_id):
    #   A -> B
    #   B -> C
    #   C -> E
    #   D -> E

    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) VALUES ('E', NULL)"  # noqa: S608
    )

    table_e_id = fetch_row_id_by_uid(
        self_referential_pg_table, ["id"], {"uid": "E"}, cursor
    )[0]

    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) VALUES ('C', {table_e_id})"  # noqa: S608
    )
    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) VALUES ('D', {table_e_id})"  # noqa: S608
    )

    table_c_id = fetch_row_id_by_uid(
        self_referential_pg_table, ["id"], {"uid": "C"}, cursor
    )[0]

    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) VALUES ('B', {table_c_id})"  # noqa: S608
    )

    table_b_id = fetch_row_id_by_uid(
        self_referential_pg_table, ["id"], {"uid": "B"}, cursor
    )[0]

    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) VALUES ('A', {table_b_id})"  # noqa: S608
    )

    _update_write_log_triggers(table_list={self_referential_pg_table})

    write_log = WriteLog.objects.all().order_by("id")

    assert len(write_log) == 5
    assert all(h.table.table_name == self_referential_pg_table for h in write_log)

    # Backfill order should be E, (D and C in some order), B, A

    assert write_log[0].record_uid == ["E"]
    assert write_log[0].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["parent_id"], {"uid": None}]
        ],
    }
    assert write_log[0].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[1].record_uid == ["D"] or write_log[1].record_uid == ["C"]
    assert write_log[1].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["parent_id"], {"uid": "E"}]
        ],
    }
    assert write_log[1].event_type == WriteLogEventTypes.INSERT.value

    assert (
        write_log[2].record_uid == ["D"] or write_log[2].record_uid == ["C"]
    ) and write_log[2].record_uid != write_log[1].record_uid
    assert write_log[2].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["parent_id"], {"uid": "E"}]
        ],
    }
    assert write_log[2].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[3].record_uid == ["B"]
    assert write_log[3].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["parent_id"], {"uid": "C"}]
        ],
    }
    assert write_log[3].event_type == WriteLogEventTypes.INSERT.value

    assert write_log[4].record_uid == ["A"]
    assert write_log[4].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [write_log[0].table.id, ["parent_id"], {"uid": "B"}]
        ],
    }
    assert write_log[4].event_type == WriteLogEventTypes.INSERT.value


@pytest.mark.pg_integration
def test_write_log_records_space_uids_properly(table_with_space_ref, fake_space):
    cursor = django_connection.cursor()

    _update_write_log_triggers(
        table_list={table_with_space_ref, "lamindb_space"},
        tables_with_installed_triggers={"lamindb_space"},
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

    write_log = WriteLog.objects.all().order_by("id")

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


@pytest.fixture(scope="function")
def drop_all_write_log_triggers_after_test():
    yield

    cursor = django_connection.cursor()

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


@pytest.mark.pg_integration
def test_write_log_install_triggers_on_existing_lamindb_models(
    drop_all_write_log_triggers_after_test,
):
    installer = PostgresWriteLogRecordingTriggerInstaller(
        connection=django_connection, db_metadata=PostgresDatabaseMetadataWrapper()
    )
    installer.update_write_log_triggers()


@pytest.fixture(scope="function")
def fake_run():
    transform = Transform("my test transform").save()
    run = Run(transform).save()

    yield run

    run.delete()


@pytest.fixture(scope="function")
def aux_artifact():
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file = Path(tmpdirname) / "foo.txt"

        with test_file.open("w") as fp:
            fp.write("hello world")

        artifact = Artifact(
            data=str(test_file),
            description="a fake aux artifact",
            kind="__lamindb_run__",  # type: ignore
            run=None,
        ).save()

        yield artifact

        artifact.delete(permanent=True)


@pytest.fixture(scope="function")
def normal_artifact(fake_run):
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file = Path(tmpdirname) / "bar.txt"

        with test_file.open("w") as fp:
            fp.write("another test file")

        artifact = Artifact(
            data=str(test_file),
            description="a fake normal artifact",
            kind="dataset",
            run=fake_run,
        ).save()

        yield artifact

        artifact.delete(permanent=True)


@pytest.mark.pg_integration
def test_aux_artifacts_backfill_before_real_ones(
    aux_artifact, normal_artifact, fake_run, drop_all_write_log_triggers_after_test
):
    installer = PostgresWriteLogRecordingTriggerInstaller(
        connection=django_connection, db_metadata=PostgresDatabaseMetadataWrapper()
    )
    installer.update_write_log_triggers()

    write_log = [
        w
        for w in WriteLog.objects.all().order_by("id")
        if w.table.table_name
        in ("lamindb_artifact", "lamindb_transform", "lamindb_run")
    ]

    assert len(write_log) == 4

    assert write_log[0].table.table_name == "lamindb_artifact"
    assert write_log[0].record_uid == [aux_artifact.uid]

    assert write_log[1].table.table_name == "lamindb_transform"
    assert write_log[1].record_uid == [normal_artifact.transform.uid]

    assert write_log[2].table.table_name == "lamindb_run"
    assert write_log[2].record_uid == [fake_run.uid]

    assert write_log[3].table.table_name == "lamindb_artifact"
    assert write_log[3].record_uid == [normal_artifact.uid]
