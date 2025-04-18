from typing import Any, Generator, Optional
from unittest import mock
from unittest.mock import MagicMock, Mock

import pytest
from django.db import connection as django_connection
from django.db import transaction
from lamindb.history._trigger_installer import (
    FOREIGN_KEYS_LIST_COLUMN_NAME,
    HistoryEventTypes,
    PostgresHistoryRecordingTriggerInstaller,
)
from lamindb.models.history import History, HistoryMigrationState, HistoryTableState


def create_spy(obj, method_name):
    original_method = getattr(obj, method_name)
    mock = Mock(wraps=original_method)
    setattr(obj, method_name, mock)
    return mock


class FakeCursor:
    def __init__(self):
        self.reset()

    def reset(self):
        self._last_result = None
        self._tables_with_triggers = []

    def _set_tables_with_triggers(self, tables: list[str]):
        self._tables_with_triggers = tables

    def execute(self, query, parameters: Optional[list[Any]] = None):
        if "DISTINCT event_object_table" in query:
            self._last_result = [(t,) for t in self._tables_with_triggers]

    def fetchall(self):
        return self._last_result


@pytest.fixture(scope="function")
def fake_cursor():
    cursor = FakeCursor()
    yield cursor
    cursor.reset()


@pytest.fixture(scope="function")
def fake_db(fake_cursor):
    db = MagicMock()
    db.cursor.return_value = fake_cursor

    return db


@pytest.fixture(scope="function", autouse=True)
def history_state():
    # Clean up after any tests that modify history
    History.objects.all().delete()
    HistoryTableState.objects.all().delete()
    HistoryMigrationState.objects.all().delete()

    yield

    History.objects.all().delete()
    HistoryTableState.objects.all().delete()
    HistoryMigrationState.objects.all().delete()


def test_get_tables_with_installed_triggers(fake_db, fake_cursor):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)
    execute = create_spy(fake_cursor, "execute")

    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b", "table_c"})
    fake_cursor._set_tables_with_triggers(["table_a", "table_c"])

    tables_with_triggers = installer.get_tables_with_installed_triggers(
        cursor=fake_cursor
    )
    execute.assert_called_once()

    assert tables_with_triggers == {
        "table_a",
        "table_c",
    }

    assert "DISTINCT event_object_table" in execute.call_args[0][0]


def test_updating_history_triggers_installs_table_state(fake_db, fake_cursor):
    assert len(set(HistoryTableState.objects.all())) == 0

    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)
    installer.install_triggers = MagicMock(return_value=None)

    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b", "table_c"})
    fake_cursor._set_tables_with_triggers(["table_a", "table_c"])

    installer.update_history_triggers()

    new_table_states = set(HistoryTableState.objects.all())

    assert len(new_table_states) == 3
    assert {ts.table_name for ts in new_table_states} == {
        "table_a",
        "table_b",
        "table_c",
    }
    assert all(ts.backfilled is False for ts in new_table_states)


def test_update_history_triggers_only_install_table_state_once(fake_db, fake_cursor):
    assert len(set(HistoryTableState.objects.all())) == 0

    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)
    installer.install_triggers = MagicMock(return_value=None)

    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b", "table_c"})
    fake_cursor._set_tables_with_triggers(["table_a", "table_c"])

    # Run update_history_triggers() twice.
    installer.update_history_triggers(update_all=True)
    installer.update_history_triggers(update_all=True)

    new_table_states = set(HistoryTableState.objects.all())

    # There should only be one table state per table, even though we've run more than once.
    assert len(new_table_states) == 3
    assert {ts.table_name for ts in new_table_states} == {
        "table_a",
        "table_b",
        "table_c",
    }
    assert all(ts.backfilled is False for ts in new_table_states)


@pytest.mark.pg_integration
def test_update_triggers_installs_migration_state():
    assert len(set(HistoryMigrationState.objects.all())) == 0

    installer = PostgresHistoryRecordingTriggerInstaller(connection=django_connection)

    installer.install_triggers = MagicMock(return_value=None)

    installer.update_history_triggers()

    new_migration_states = set(HistoryMigrationState.objects.all())

    assert len(new_migration_states) == 1

    # Updating triggers a second time shouldn't add more migration state, since
    # no migrations have been performed between the two updates.

    installer.update_history_triggers()

    new_migration_states = set(HistoryMigrationState.objects.all())

    assert len(new_migration_states) == 1


def test_update_history_triggers_skips_existing_triggers(fake_db, fake_cursor):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)

    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b", "table_c"})
    fake_cursor._set_tables_with_triggers(["table_a", "table_c"])

    installer.install_triggers = MagicMock(return_value=None)

    installer.update_history_triggers()

    installer._get_db_tables.assert_called_once()
    installer.install_triggers.assert_called_once_with("table_b", fake_cursor)


@pytest.fixture(scope="function")
def no_uid_pg_table(history_state) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS history_table_test_no_uid
(id SERIAL PRIMARY KEY, str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield "history_table_test_no_uid"

    cursor.execute("DROP TABLE IF EXISTS history_table_test_no_uid")


@pytest.mark.pg_integration
def test_foreign_key_to_table_without_uid_fails(no_uid_pg_table):
    with pytest.raises(ValueError):
        _update_history_triggers({no_uid_pg_table})


def test_sql_injectable_table_names_fail(fake_cursor):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)
    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b"})

    with pytest.raises(ValueError):
        installer.install_triggers("bad_tabl; DROP TABLE foo", fake_cursor)


@pytest.fixture(scope="function")
def simple_pg_table(history_state) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS history_table_test_a
(id SERIAL PRIMARY KEY, uid VARCHAR(8), str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield "history_table_test_a"

    cursor.execute("DROP TABLE IF EXISTS history_table_test_a")


def _update_history_triggers(table_list: set[str]):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=django_connection)

    with mock.patch.object(installer, "_get_db_tables") as get_db_tables:
        get_db_tables.return_value = table_list

        installer.update_history_triggers()


@pytest.mark.pg_integration
def test_simple_table_trigger(simple_pg_table: str):
    _update_history_triggers({simple_pg_table})

    cursor = django_connection.cursor()

    cursor.execute("SELECT * FROM pg_trigger WHERE tgname LIKE 'lamindb_history_%';")

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

    history = History.objects.all().order_by("seqno")

    assert len(history) == 4

    assert [h.table.table_name == simple_pg_table for h in history]

    assert history[0].record_uid == "abc123"
    assert history[0].record_data == {
        "uid": "abc123",
        "int_col": 42,
        "str_col": "hello world",
        "bool_col": False,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert history[0].event_type == HistoryEventTypes.INSERT.value

    assert history[1].record_uid == "def456"
    assert history[1].record_data == {
        "uid": "def456",
        "int_col": 99,
        "str_col": "another row",
        "bool_col": True,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert history[1].event_type == HistoryEventTypes.INSERT.value

    assert history[2].record_uid == "def456"
    assert history[2].record_data == {
        "uid": "def456",
        "int_col": 22,
        "str_col": "another row",
        "bool_col": True,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert history[2].event_type == HistoryEventTypes.UPDATE.value

    assert history[3].record_uid == "abc123"
    assert history[3].record_data is None
    assert history[3].event_type == HistoryEventTypes.DELETE.value


@pytest.fixture(scope="function")
def foreignkey_pg_table(simple_pg_table) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS history_table_test_b
(
    id SERIAL PRIMARY KEY, uid VARCHAR(8), table_a_id int,
    CONSTRAINT fk_table_a FOREIGN KEY(table_a_id) REFERENCES {simple_pg_table}(id)
    ON DELETE CASCADE
);
""")

    yield "history_table_test_b"

    cursor.execute("DROP TABLE IF EXISTS history_table_test_b")


@pytest.mark.pg_integration
def test_triggers_with_foreign_keys(simple_pg_table, foreignkey_pg_table):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=django_connection)
    installer._get_db_tables = MagicMock(
        return_value={simple_pg_table, foreignkey_pg_table}
    )

    installer.update_history_triggers()

    cursor = django_connection.cursor()

    with transaction.atomic():
        cursor.execute(
            f"INSERT INTO {simple_pg_table} (uid, str_col, bool_col, int_col) "  # noqa: S608
            "VALUES ('abc123', 'hello world', false, 42)"
        )
        cursor.execute(f"SELECT id FROM {simple_pg_table} WHERE uid = 'abc123'")  # noqa: S608
        table_a_row_id = cursor.fetchone()[0]

        cursor.execute(
            f"INSERT INTO {foreignkey_pg_table} (uid, table_a_id) VALUES ('foo333', {table_a_row_id})"  # noqa: S608
        )
    cursor.execute(
        f"INSERT INTO {foreignkey_pg_table} (uid, table_a_id) VALUES ('bar444', NULL)"  # noqa: S608
    )
    cursor.execute(f"DELETE FROM {simple_pg_table} WHERE uid = 'abc123'")  # noqa: S608

    history = History.objects.all().order_by("seqno")

    assert len(history) == 5

    assert history[0].table.table_name == simple_pg_table
    assert history[0].record_uid == "abc123"
    assert history[0].record_data == {
        "uid": "abc123",
        "int_col": 42,
        "str_col": "hello world",
        "bool_col": False,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert history[0].event_type == HistoryEventTypes.INSERT.value

    assert history[1].table.table_name == foreignkey_pg_table
    assert history[1].record_uid == "foo333"
    assert history[1].record_data == {
        "uid": "foo333",
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [history[0].table.id, ["table_a_id"], {"uid": "abc123"}]
        ],
    }
    assert history[1].event_type == HistoryEventTypes.INSERT.value

    assert history[2].table.table_name == foreignkey_pg_table
    assert history[2].record_uid == "bar444"
    assert history[2].event_type == HistoryEventTypes.INSERT.value
    assert history[2].record_data == {
        "uid": "bar444",
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [history[0].table_id, ["table_a_id"], {"uid": None}]
        ],
    }

    assert history[3].table.table_name == simple_pg_table
    assert history[3].record_uid == "abc123"
    assert history[3].event_type == HistoryEventTypes.DELETE.value
    assert history[3].record_data is None

    assert history[4].table.table_name == foreignkey_pg_table
    assert history[4].record_uid == "foo333"
    assert history[4].event_type == HistoryEventTypes.DELETE.value
    assert history[4].record_data is None


@pytest.fixture(scope="function")
def self_referential_pg_table() -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS history_table_test_c
(
    id SERIAL PRIMARY KEY, uid VARCHAR(8), parent_id int,
    CONSTRAINT fk_parent FOREIGN KEY(parent_id) REFERENCES history_table_test_c(id)
    ON DELETE CASCADE
);
""")

    yield "history_table_test_c"

    cursor.execute("DROP TABLE IF EXISTS history_table_test_c")


@pytest.mark.pg_integration
def test_triggers_with_self_references(self_referential_pg_table):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=django_connection)
    installer._get_db_tables = MagicMock(return_value={self_referential_pg_table})

    installer.update_history_triggers()

    cursor = django_connection.cursor()

    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) "  # noqa: S608
        "VALUES ('abc123', NULL)"
    )
    cursor.execute(f"SELECT id FROM {self_referential_pg_table} WHERE uid = 'abc123'")  # noqa: S608
    parent_row_id = cursor.fetchone()[0]
    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) VALUES ('def345', {parent_row_id})"  # noqa: S608
    )

    history = History.objects.all().order_by("seqno")

    assert len(history) == 2

    assert history[0].table.table_name == self_referential_pg_table
    assert history[0].record_uid == "abc123"
    assert history[0].record_data == {
        "uid": "abc123",
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [history[0].table.id, ["parent_id"], {"uid": None}]
        ],
    }
    assert history[0].event_type == HistoryEventTypes.INSERT.value

    assert history[1].table.table_name == self_referential_pg_table
    assert history[1].record_uid == "def345"
    assert history[1].record_data == {
        "uid": "def345",
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [history[0].table.id, ["parent_id"], {"uid": "abc123"}]
        ],
    }
    assert history[1].event_type == HistoryEventTypes.INSERT.value


@pytest.fixture(scope="function")
def composite_primary_key_table() -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS history_table_test_d
(
    key_a INT,
    key_b INT,
    uid VARCHAR(8),
    PRIMARY KEY (key_a, key_b)
);
""")

    yield "history_table_test_d"

    cursor.execute("DROP TABLE IF EXISTS history_table_test_d")


@pytest.fixture(scope="function")
def composite_foreign_key_table(
    composite_primary_key_table,
) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS history_table_test_e
(
    id SERIAL PRIMARY KEY,
    table_d_key_a INT,
    table_d_key_b INT,
    uid VARCHAR(8),
    CONSTRAINT fk_parent FOREIGN KEY(table_d_key_a, table_d_key_b) REFERENCES {composite_primary_key_table} (key_a, key_b)
);
""")

    yield "history_table_test_e"

    cursor.execute("DROP TABLE IF EXISTS history_table_test_e")


@pytest.mark.pg_integration
def test_triggers_with_composite_primary_key(
    composite_primary_key_table, composite_foreign_key_table
):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=django_connection)
    installer._get_db_tables = MagicMock(
        return_value={composite_primary_key_table, composite_foreign_key_table}
    )

    installer.update_history_triggers()

    cursor = django_connection.cursor()

    cursor.execute(
        f"INSERT INTO {composite_primary_key_table} (key_a, key_b, uid) "  # noqa: S608
        "VALUES (42, 21, 'xyz123')"
    )
    cursor.execute(
        f"INSERT INTO {composite_foreign_key_table} (table_d_key_a, table_d_key_b, uid) "  # noqa: S608
        f"VALUES (42, 21, 'qrs234')"
    )

    history = History.objects.all().order_by("seqno")

    assert len(history) == 2

    assert history[0].table.table_name == composite_primary_key_table
    assert history[0].record_uid == "xyz123"
    assert history[0].record_data == {
        "uid": "xyz123",
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert history[0].event_type == HistoryEventTypes.INSERT.value

    assert history[1].table.table_name == composite_foreign_key_table
    assert history[1].record_uid == "qrs234"
    assert history[1].record_data == {
        "uid": "qrs234",
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [history[0].table.id, ["table_d_key_a", "table_d_key_b"], {"uid": "xyz123"}]
        ],
    }
    assert history[1].event_type == HistoryEventTypes.INSERT.value
