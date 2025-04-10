from unittest.mock import MagicMock, Mock

import pytest
from lamindb.core.history._trigger_installer import (
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

    def execute(self, query, parameters):
        if "DISTINCT event_object_table" in query:
            self._last_result = [("table_a",), ("table_b",), ("table_c",)]

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
    db.cursor.returnvalue = fake_cursor

    return db


def test_get_tables_with_installed_triggers(fake_db, fake_cursor):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)
    execute = create_spy(fake_cursor, "execute")

    execute.assert_called_once()

    assert installer.get_tables_with_installed_triggers() == {
        "table_a",
        "table_b",
        "table_c",
    }

    assert "DISTINCT event_object_table" in execute.call_args[0]
