from unittest.mock import MagicMock, Mock

import pytest
from lamindb.history._trigger_installer import (
    PostgresHistoryRecordingTriggerInstaller,
)


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

    tables_with_triggers = installer.get_tables_with_installed_triggers()
    execute.assert_called_once()

    assert tables_with_triggers == {
        "table_a",
        "table_c",
    }

    assert "DISTINCT event_object_table" in execute.call_args[0][0]


def test_update_history_triggers_skips_existing_triggers(fake_db, fake_cursor):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)

    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b", "table_c"})
    installer.install_triggers = MagicMock(return_value=None)

    installer.update_history_triggers()

    installer._get_db_tables.assert_called_once()
    installer.install_triggers.assert_called_once_with("table_b")
