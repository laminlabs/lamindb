import json
from typing import TYPE_CHECKING, Generator, cast
from unittest.mock import ANY, MagicMock

import pytest
from django.db import connection as django_connection_proxy
from django.db import transaction
from django.db.backends.utils import CursorWrapper
from lamindb.history._db_metadata_wrapper import (
    PostgresDatabaseMetadataWrapper,
)
from lamindb.history._trigger_installer import (
    FOREIGN_KEYS_LIST_COLUMN_NAME,
    HistoryEventTypes,
    PostgresHistoryRecordingTriggerInstaller,
)
from lamindb.models.history import History, HistoryMigrationState, HistoryTableState
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


@pytest.fixture(scope="function")
def table_a(history_state) -> Generator[str, None, None]:
    table_name = "history_test_table_a"
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS {table_name}
(id SERIAL PRIMARY KEY, str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield table_name

    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")


@pytest.fixture(scope="function")
def table_b(history_state) -> Generator[str, None, None]:
    table_name = "history_test_table_b"
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS {table_name}
(id SERIAL PRIMARY KEY, str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield table_name

    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")


@pytest.fixture(scope="function")
def table_c(history_state) -> Generator[str, None, None]:
    table_name = "history_test_table_c"
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS {table_name}
(id SERIAL PRIMARY KEY, str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield table_name

    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")


@pytest.mark.pg_integration
def test_updating_history_triggers_installs_table_state(table_a, table_b, table_c):
    assert len(set(HistoryTableState.objects.all())) == 0

    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b, table_c})
    fake_db_metadata.set_tables_with_installed_triggers({table_a, table_c})

    installer = PostgresHistoryRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )
    installer.install_triggers = MagicMock(return_value=None)

    installer.update_history_triggers()

    new_table_states = set(HistoryTableState.objects.all())

    assert len(new_table_states) == 3
    assert {ts.table_name for ts in new_table_states} == {
        table_a,
        table_b,
        table_c,
    }
    assert all(ts.backfilled is False for ts in new_table_states)


@pytest.mark.pg_integration
def test_update_history_triggers_only_install_table_state_once(
    table_a, table_b, table_c
):
    assert len(set(HistoryTableState.objects.all())) == 0

    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b, table_c})
    fake_db_metadata.set_tables_with_installed_triggers({table_a, table_c})

    installer = PostgresHistoryRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )
    installer.install_triggers = MagicMock(return_value=None)

    # Run update_history_triggers() twice.
    installer.update_history_triggers(update_all=True)
    installer.update_history_triggers(update_all=True)

    new_table_states = set(HistoryTableState.objects.all())

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
    assert len(set(HistoryMigrationState.objects.all())) == 0

    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b, table_c})
    fake_db_metadata.set_tables_with_installed_triggers(set())

    installer = PostgresHistoryRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )

    installer.install_triggers = MagicMock(return_value=None)

    installer.update_history_triggers()

    new_migration_states = set(HistoryMigrationState.objects.all())

    assert len(new_migration_states) == 1

    # Updating triggers a second time shouldn't add more migration state, since
    # no migrations have been performed between the two updates.

    installer.update_history_triggers()

    new_migration_states = set(HistoryMigrationState.objects.all())

    assert len(new_migration_states) == 1


@pytest.mark.pg_integration
def test_update_history_triggers_skips_existing_triggers(table_a, table_b, table_c):
    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b, table_c})
    fake_db_metadata.set_tables_with_installed_triggers({table_a, table_c})

    installer = PostgresHistoryRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )

    installer.install_triggers = MagicMock(return_value=None)

    installer.update_history_triggers()

    installer.install_triggers.assert_called_once_with(table_b, ANY)


@pytest.fixture(scope="function")
def no_uid_pg_table(history_state) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE IF NOT EXISTS history_table_test_no_uid
(id SERIAL PRIMARY KEY, str_col VARCHAR(50), bool_col BOOLEAN, int_col INT);
""")

    yield "history_table_test_no_uid"

    cursor.execute("DROP TABLE IF EXISTS history_table_test_no_uid")


@pytest.fixture(scope="function")
def foreign_key_to_no_uid_table(
    history_state, no_uid_pg_table
) -> Generator[str, None, None]:
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE IF NOT EXISTS history_table_test_no_uid_fk
(id SERIAL PRIMARY KEY, primary_id INT, CONSTRAINT fk FOREIGN KEY (primary_id) REFERENCES {no_uid_pg_table}(id));
""")

    yield "history_table_test_no_uid_fk"

    cursor.execute("DROP TABLE IF EXISTS history_table_test_no_uid_fk")


@pytest.mark.pg_integration
def test_foreign_key_to_table_without_uid_fails(
    no_uid_pg_table, foreign_key_to_no_uid_table
):
    with pytest.raises(ValueError):
        _update_history_triggers({no_uid_pg_table, foreign_key_to_no_uid_table})


@pytest.mark.pg_integration
def test_sql_injectable_table_names_fail(table_a, table_b):
    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables({table_a, table_b})
    fake_db_metadata.set_tables_with_installed_triggers(set())

    installer = PostgresHistoryRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )

    with pytest.raises(ValueError):
        installer.install_triggers(
            "bad_tabl; DROP TABLE foo", django_connection.cursor()
        )


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
    fake_db_metadata = FakeMetadataWrapper()
    fake_db_metadata.set_db_tables(table_list)
    fake_db_metadata.set_tables_with_installed_triggers(set())

    installer = PostgresHistoryRecordingTriggerInstaller(
        connection=django_connection, db_metadata=fake_db_metadata
    )

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

    assert history[0].record_uid == json.dumps(["abc123"])
    assert history[0].record_data == {
        "int_col": 42,
        "str_col": "hello world",
        "bool_col": False,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert history[0].event_type == HistoryEventTypes.INSERT.value

    assert history[1].record_uid == json.dumps(["def456"])
    assert history[1].record_data == {
        "int_col": 99,
        "str_col": "another row",
        "bool_col": True,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert history[1].event_type == HistoryEventTypes.INSERT.value

    assert history[2].record_uid == json.dumps(["def456"])
    assert history[2].record_data == {
        "int_col": 22,
        "str_col": "another row",
        "bool_col": True,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert history[2].event_type == HistoryEventTypes.UPDATE.value

    assert history[3].record_uid == json.dumps(["abc123"])
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
    _update_history_triggers({simple_pg_table, foreignkey_pg_table})

    cursor = django_connection.cursor()

    with transaction.atomic():
        cursor.execute(
            f"INSERT INTO {simple_pg_table} (uid, str_col, bool_col, int_col) "  # noqa: S608
            "VALUES ('abc123', 'hello world', false, 42)"
        )
        cursor.execute(f"SELECT id FROM {simple_pg_table} WHERE uid = 'abc123'")  # noqa: S608
        table_a_row = cursor.fetchone()
        assert table_a_row is not None
        table_a_row_id = table_a_row[0]

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
    assert history[0].record_uid == json.dumps(["abc123"])
    assert history[0].record_data == {
        "int_col": 42,
        "str_col": "hello world",
        "bool_col": False,
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert history[0].event_type == HistoryEventTypes.INSERT.value

    assert history[1].table.table_name == foreignkey_pg_table
    assert history[1].record_uid == json.dumps(["foo333"])
    assert history[1].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [history[0].table.id, ["table_a_id"], {"uid": "abc123"}]
        ],
    }
    assert history[1].event_type == HistoryEventTypes.INSERT.value

    assert history[2].table.table_name == foreignkey_pg_table
    assert history[2].record_uid == json.dumps(["bar444"])
    assert history[2].event_type == HistoryEventTypes.INSERT.value
    assert history[2].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [history[0].table.id, ["table_a_id"], {"uid": None}]
        ],
    }

    assert history[3].table.table_name == simple_pg_table
    assert history[3].record_uid == json.dumps(["abc123"])
    assert history[3].event_type == HistoryEventTypes.DELETE.value
    assert history[3].record_data is None

    assert history[4].table.table_name == foreignkey_pg_table
    assert history[4].record_uid == json.dumps(["foo333"])
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
    _update_history_triggers({self_referential_pg_table})

    cursor = django_connection.cursor()

    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) "  # noqa: S608
        "VALUES ('abc123', NULL)"
    )
    cursor.execute(f"SELECT id FROM {self_referential_pg_table} WHERE uid = 'abc123'")  # noqa: S608

    parent_row = cursor.fetchone()
    assert parent_row is not None
    parent_row_id = parent_row[0]

    cursor.execute(
        f"INSERT INTO {self_referential_pg_table} (uid, parent_id) VALUES ('def345', {parent_row_id})"  # noqa: S608
    )

    history = History.objects.all().order_by("seqno")

    assert len(history) == 2

    assert history[0].table.table_name == self_referential_pg_table
    assert history[0].record_uid == json.dumps(["abc123"])
    assert history[0].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [history[0].table.id, ["parent_id"], {"uid": None}]
        ],
    }
    assert history[0].event_type == HistoryEventTypes.INSERT.value

    assert history[1].table.table_name == self_referential_pg_table
    assert history[1].record_uid == json.dumps(["def345"])
    assert history[1].record_data == {
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
    _update_history_triggers({composite_primary_key_table, composite_foreign_key_table})

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
    assert history[0].record_uid == json.dumps(["xyz123"])
    assert history[0].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [],
    }
    assert history[0].event_type == HistoryEventTypes.INSERT.value

    assert history[1].table.table_name == composite_foreign_key_table
    assert history[1].record_uid == json.dumps(["qrs234"])
    assert history[1].record_data == {
        FOREIGN_KEYS_LIST_COLUMN_NAME: [
            [history[0].table.id, ["table_d_key_a", "table_d_key_b"], {"uid": "xyz123"}]
        ],
    }
    assert history[1].event_type == HistoryEventTypes.INSERT.value
