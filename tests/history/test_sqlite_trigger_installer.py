from typing import Any, Optional
from unittest.mock import MagicMock

import pytest
from lamindb.history._trigger_installer import SQLiteHistoryRecordingTriggerInstaller


class FakeCursor:
    def __init__(self):
        self.reset()

    def execute(self, query, parameters: Optional[list[Any]] = None):
        pass

    def fetchall(self):
        return self._last_result

    def reset(self):
        self._last_result = None
        self._tables_with_triggers = []
        self._constraints = []
        self._column_names = {}


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
    # TODO: replace this with a real test once the installer is implemented

    installer = SQLiteHistoryRecordingTriggerInstaller(fake_db)

    assert installer.get_tables_with_installed_triggers(fake_cursor) == set()


def test_install_triggers(fake_db, fake_cursor):
    # TODO: replace this with a real test once the installer is implemented

    installer = SQLiteHistoryRecordingTriggerInstaller(fake_db)

    installer.install_triggers("table_a", fake_cursor)
