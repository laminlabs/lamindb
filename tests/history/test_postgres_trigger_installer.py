from unittest.mock import MagicMock, Mock

import pytest
from django.db import connection as django_connection
from lamindb.history._trigger_installer import (
    PostgresHistoryRecordingTriggerInstaller,
)
from lamindb.models.history import HistoryMigrationState, HistoryTableState


def create_spy(obj, method_name):
    original_method = getattr(obj, method_name)
    mock = Mock(wraps=original_method)
    setattr(obj, method_name, mock)
    return mock


class FakeCursor:
    def __init__(self):
        self._last_result = None

    def execute(self, query, parameters=None):
        if "DISTINCT event_object_table" in query:
            self._last_result = [("table_a",), ("table_c",)]

    def fetchall(self):
        return self._last_result

    def reset(self):
        self._last_result = None


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


def test_get_tables_with_installed_triggers(fake_db, fake_cursor):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)
    execute = create_spy(fake_cursor, "execute")

    tables_with_triggers = installer.get_tables_with_installed_triggers(
        cursor=fake_cursor
    )
    execute.assert_called_once()

    assert tables_with_triggers == {
        "table_a",
        "table_c",
    }

    assert "DISTINCT event_object_table" in execute.call_args[0][0]


@pytest.fixture(scope="function", autouse=True)
def history_table_state():
    # Make sure we're not polluting other tests with HistoryTableState created by
    # any given test.
    HistoryTableState.objects.all().delete()

    yield

    HistoryTableState.objects.all().delete()


@pytest.fixture(scope="function", autouse=True)
def history_migration_state():
    # Make sure we're not polluting other tests with HistoryMigrationState created by
    # any given test.
    HistoryMigrationState.objects.all().delete()

    yield set()

    HistoryMigrationState.objects.all().delete()


def test_updating_history_triggers_installs_table_state(fake_db):
    assert len(set(HistoryTableState.objects.all())) == 0

    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)

    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b", "table_c"})
    installer.install_triggers = MagicMock(return_value=None)

    installer.update_history_triggers()

    new_table_states = set(HistoryTableState.objects.all())

    assert len(new_table_states) == 3
    assert {ts.table_name for ts in new_table_states} == {
        "table_a",
        "table_b",
        "table_c",
    }
    assert all(ts.backfilled is False for ts in new_table_states)


def test_update_history_triggers_only_install_table_state_once(fake_db):
    assert len(set(HistoryTableState.objects.all())) == 0

    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)

    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b", "table_c"})
    installer.install_triggers = MagicMock(return_value=None)

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


def test_update_triggers_installs_migration_state():
    assert len(set(HistoryMigrationState.objects.all())) == 0

    installer = PostgresHistoryRecordingTriggerInstaller(connection=django_connection)

    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b", "table_c"})
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
    installer.install_triggers = MagicMock(return_value=None)

    installer.update_history_triggers()

    installer._get_db_tables.assert_called_once()
    installer.install_triggers.assert_called_once_with("table_b", fake_cursor)
