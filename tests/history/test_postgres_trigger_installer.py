import re
from collections import namedtuple
from typing import Any, Literal, Optional
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


Constraint = namedtuple(
    "Constraint", "source_table_name type column_name target_table_name"
)


class FakeCursor:
    def __init__(self):
        self.reset()

    def _add_constraint(
        self,
        table_name,
        constraint_type: Literal["PRIMARY KEY", "FOREIGN KEY"],
        column_name: str,
        target_table_name: str,
    ):
        self._constraints.append(
            Constraint(table_name, constraint_type, column_name, target_table_name)
        )

    def _set_column_names(self, table_name: str, column_names: list[str]):
        self._column_names[table_name] = column_names

    def _set_tables_with_triggers(self, tables: list[str]):
        self._tables_with_triggers = tables

    def execute(self, query, parameters: Optional[list[Any]] = None):
        if "DISTINCT event_object_table" in query:
            self._last_result = [(t,) for t in self._tables_with_triggers]
        elif "tc.table_name = %s" in query:
            self._last_result = [
                c[1:] for c in self._constraints if c.source_table_name == parameters[0]
            ]
        elif "SELECT column_name FROM information_schema.columns" in query:
            self._last_result = [(c,) for c in self._column_names[parameters[0]]]

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


def test_install_triggers_no_foreign_keys(fake_db, fake_cursor):
    execute = create_spy(fake_cursor, "execute")

    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)
    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b", "table_c"})

    fake_cursor._set_tables_with_triggers(["table_b", "table_c"])
    fake_cursor._add_constraint("table_a", "PRIMARY KEY", "id", "table_a")
    fake_cursor._set_column_names("table_a", ["id", "uid", "foo", "bar", "baz"])

    installer.update_history_triggers()

    create_trigger_calls = [
        c for c in execute.call_args_list if "CREATE OR REPLACE TRIGGER" in c.args[0]
    ]

    assert len(create_trigger_calls) == 3

    assert len([c for c in create_trigger_calls if "AFTER INSERT" in c.args[0]]) == 1
    assert len([c for c in create_trigger_calls if "AFTER UPDATE" in c.args[0]]) == 1
    assert len([c for c in create_trigger_calls if "AFTER DELETE" in c.args[0]]) == 1

    # Ideally we'd parse the trigger with something like SQLGlot but this'll do for now.

    create_function_calls = [
        c for c in execute.call_args_list if "CREATE OR REPLACE FUNCTION" in c.args[0]
    ]

    assert len(create_function_calls) == 1

    create_function_sql = re.sub(
        r"\s+", " ", create_function_calls[0].args[0].replace("\n", "")
    )

    declarations = [
        l for l in create_function_sql.split(";") if l.startswith("DECLARE")
    ]

    # Make sure we're referencing the right table when declaring table_id
    assert any("table_id smallint :=" in d and "'table_a'" in d for d in declarations)


def test_install_triggers_with_foreign_keys(fake_db, fake_cursor):
    execute = create_spy(fake_cursor, "execute")

    installer = PostgresHistoryRecordingTriggerInstaller(connection=fake_db)
    installer._get_db_tables = MagicMock(return_value={"table_a", "table_b"})

    fake_cursor._add_constraint("table_a", "PRIMARY KEY", "id", "table_a")
    fake_cursor._add_constraint("table_a", "FOREIGN KEY", "table_b_id", "table_b")
    fake_cursor._set_column_names("table_a", ["id", "uid", "table_b_id", "foo"])
    fake_cursor._set_column_names("table_b", ["id", "uid", "bar"])
    fake_cursor._set_tables_with_triggers(["table_b"])

    installer.update_history_triggers()

    create_trigger_calls = [
        c for c in execute.call_args_list if "CREATE OR REPLACE TRIGGER" in c.args[0]
    ]

    assert len(create_trigger_calls) == 3

    assert len([c for c in create_trigger_calls if "AFTER INSERT" in c.args[0]]) == 1
    assert len([c for c in create_trigger_calls if "AFTER UPDATE" in c.args[0]]) == 1
    assert len([c for c in create_trigger_calls if "AFTER DELETE" in c.args[0]]) == 1

    create_function_calls = [
        c for c in execute.call_args_list if "CREATE OR REPLACE FUNCTION" in c.args[0]
    ]

    assert len(create_function_calls) == 1

    create_function_sql = re.sub(
        r"\s+", " ", create_function_calls[0].args[0].replace("\n", "")
    )

    declarations = [
        l for l in create_function_sql.split(";") if l.startswith("DECLARE")
    ]

    # We should be declaring a variable that's extracting the UID from table_b
    assert any(
        d.startswith("DECLARE fkey_table_b_id_uid") and "SELECT uid FROM table_b" in d
        for d in declarations
    )

    # We should be adding the declared variable to jsonb_build_object someplace, with a marker
    # on the object's key to indicate that it's a UID reference
    assert "'table_b_id._uid', fkey_table_b_id_uid" in create_function_sql
