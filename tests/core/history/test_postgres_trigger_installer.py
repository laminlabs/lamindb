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
def cursor():
    cursor = FakeCursor()
    yield cursor
    cursor.reset()


@pytest.fixture(scope="function")
def db(cursor):
    db = MagicMock()
    db.cursor.returnvalue = cursor

    return db


def test_get_tables_with_installed_triggers(db, cursor):
    installer = PostgresHistoryRecordingTriggerInstaller(connection=db)
    execute = create_spy(cursor, "execute")

    assert installer.get_tables_with_installed_triggers() == {
        "table_a",
        "table_b",
        "table_c",
    }

    assert execute.called_once()
    assert "DISTINCT event_object_table" in execute.call_args[0]
