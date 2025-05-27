import datetime
from typing import TYPE_CHECKING, cast

import pytest
from django.db import connection
from django.db import connection as django_connection_proxy
from django.db.backends.utils import CursorWrapper
from lamindb.core.writelog._constants import FOREIGN_KEYS_LIST_COLUMN_NAME
from lamindb.core.writelog._db_metadata_wrapper import SQLiteDatabaseMetadataWrapper
from lamindb.core.writelog._replayer import WriteLogReplayer
from lamindb.core.writelog._trigger_installer import WriteLogEventTypes
from lamindb.core.writelog._types import UIDColumns
from lamindb.core.writelog._utils import (
    get_latest_migration_state,
    update_migration_state,
    update_write_log_table_state,
)
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


class FakeMetadataWrapper(SQLiteDatabaseMetadataWrapper):
    """A fake DB metadata wrapper that allows us to control which database tables the installer will see and target."""

    def __init__(self):
        super().__init__()
        self._tables_with_triggers = set()
        self._db_tables: set[str] = set()
        self._many_to_many_tables: set[str] = set()
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


@pytest.fixture(scope="function")
def write_log_state():
    WriteLog.objects.all().delete()
    WriteLogTableState.objects.all().delete()
    WriteLogMigrationState.objects.all().delete()

    yield

    WriteLog.objects.all().delete()
    WriteLogTableState.objects.all().delete()
    WriteLogMigrationState.objects.all().delete()


def test_connection_is_sqlite():
    assert connection.vendor == "sqlite"


@pytest.fixture(scope="function")
def simple_table(write_log_state):
    cursor = django_connection.cursor()

    cursor.execute("""
CREATE TABLE simple_table (
    id integer NOT NULL PRIMARY KEY AUTOINCREMENT,
    uid varchar(20) NOT NULL UNIQUE,
    bool_col bool,
    text_col TEXT,
    timestamp_col datetime,
    date_col date,
    float_col REAL
)
""")

    yield "simple_table"

    cursor.execute("DROP TABLE IF EXISTS simple_table")


@pytest.fixture(scope="function")
def write_log_lock():
    write_log_lock = WriteLogLock.load()

    assert write_log_lock is not None

    write_log_lock.lock()

    yield write_log_lock

    write_log_lock.unlock()


def test_replayer_happy_path(simple_table, write_log_lock, write_log_state):
    db_metadata = FakeMetadataWrapper()
    db_metadata.set_db_tables({simple_table})

    cursor = django_connection.cursor()

    update_write_log_table_state({simple_table})
    update_migration_state()

    current_migration_state = get_latest_migration_state()

    simple_table_state = WriteLogTableState.objects.get(table_name=simple_table)

    replayer = WriteLogReplayer(db_metadata=db_metadata, cursor=cursor)

    input_write_log = [
        WriteLog(
            migration_state=current_migration_state,
            table=simple_table_state,
            uid="Hist001",
            record_uid=["SimpleRecord1"],
            record_data={
                "bool_col": True,
                "text_col": "hello world",
                "timestamp_col": "2025-05-23T02:03:16.913425+00:00",
                "date_col": "2025-05-23",
                "float_col": 8.675309,
                FOREIGN_KEYS_LIST_COLUMN_NAME: [],
            },
            event_type=WriteLogEventTypes.INSERT.value,
            created_at=datetime.datetime(2025, 5, 23, 12, 34, 56),
        ),
        WriteLog(
            migration_state=current_migration_state,
            table=simple_table_state,
            uid="Hist002",
            record_uid=["SimpleRecord2"],
            record_data={
                "bool_col": False,
                "text_col": "Hallo, Welt!",
                "timestamp_col": "2025-05-23T02:03:20.392310+00:00",
                "date_col": "2025-05-23",
                "float_col": 1.8007777777,
                FOREIGN_KEYS_LIST_COLUMN_NAME: [],
            },
            event_type=WriteLogEventTypes.INSERT.value,
            created_at=datetime.datetime(2025, 5, 23, 12, 41, 00),
        ),
        WriteLog(
            migration_state=current_migration_state,
            table=simple_table_state,
            uid="Hist003",
            record_uid=["SimpleRecord1"],
            record_data={
                "bool_col": False,
                "text_col": "hello world",
                "timestamp_col": "2025-05-23T02:03:16.913425+00:00",
                "date_col": "2025-05-24",
                "float_col": 8.675309,
                FOREIGN_KEYS_LIST_COLUMN_NAME: [],
            },
            event_type=WriteLogEventTypes.UPDATE.value,
            created_at=datetime.datetime(2025, 5, 23, 12, 55, 00),
        ),
        WriteLog(
            migration_state=current_migration_state,
            table=simple_table_state,
            uid="Hist004",
            record_uid=["SimpleRecord2"],
            record_data=None,
            event_type=WriteLogEventTypes.DELETE.value,
            created_at=datetime.datetime(2025, 5, 23, 12, 56, 00),
        ),
    ]

    for write_log_entry in input_write_log:
        replayer.replay(write_log_entry)
        write_log_entry.save()

    write_log = WriteLog.objects.all().order_by("id")

    assert len(write_log) == 4
    assert list(write_log) == input_write_log

    cursor.execute(f"SELECT uid, bool_col FROM {simple_table} ORDER BY uid ASC")  # noqa: S608

    rows = cursor.fetchall()
    assert len(rows) == 1
    assert rows[0][0] == "SimpleRecord1"
    assert rows[0][1] is False


@pytest.fixture(scope="function")
def many_to_many_table(simple_table, write_log_state):
    cursor = django_connection.cursor()

    cursor.execute(f"""
CREATE TABLE many_to_many_table (
    id integer NOT NULL PRIMARY KEY AUTOINCREMENT,
    simple_a_id integer,
    simple_b_id integer,
    FOREIGN KEY (simple_a_id) REFERENCES {simple_table}(id),
    FOREIGN KEY (simple_b_id) REFERENCES {simple_table}(id)
)
""")

    yield "many_to_many_table"

    cursor.execute("DROP TABLE IF EXISTS many_to_many_table")


def test_replayer_many_to_many(
    simple_table, many_to_many_table, write_log_lock, write_log_state
):
    db_metadata = FakeMetadataWrapper()
    db_metadata.set_db_tables({simple_table, many_to_many_table})
    db_metadata.set_many_to_many_db_tables({many_to_many_table})

    cursor = django_connection.cursor()

    update_write_log_table_state({simple_table, many_to_many_table})
    update_migration_state()

    current_migration_state = get_latest_migration_state()

    simple_table_state = WriteLogTableState.objects.get(table_name=simple_table)
    many_to_many_table_state = WriteLogTableState.objects.get(
        table_name=many_to_many_table
    )

    replayer = WriteLogReplayer(db_metadata=db_metadata, cursor=cursor)

    input_write_log = [
        WriteLog(
            migration_state=current_migration_state,
            table=simple_table_state,
            uid="Hist001",
            record_uid=["SimpleRecord1"],
            record_data={
                "bool_col": True,
                "text_col": "hello world",
                "timestamp_col": "2025-05-23T02:03:16.913425+00:00",
                "date_col": "2025-05-23",
                "float_col": 8.675309,
                FOREIGN_KEYS_LIST_COLUMN_NAME: [],
            },
            event_type=WriteLogEventTypes.INSERT.value,
            created_at=datetime.datetime(
                2025, 5, 23, 12, 34, 56, tzinfo=datetime.timezone.utc
            ),
        ),
        WriteLog(
            migration_state=current_migration_state,
            table=simple_table_state,
            uid="Hist002",
            record_uid=["SimpleRecord2"],
            record_data={
                "bool_col": False,
                "text_col": "Hallo, Welt!",
                "timestamp_col": "2025-05-23T02:03:20.392310+00:00",
                "date_col": "2025-05-23",
                "float_col": 1.8007777777,
                FOREIGN_KEYS_LIST_COLUMN_NAME: [],
            },
            event_type=WriteLogEventTypes.INSERT.value,
            created_at=datetime.datetime(
                2025, 5, 23, 12, 41, 00, tzinfo=datetime.timezone.utc
            ),
        ),
        WriteLog(
            migration_state=current_migration_state,
            table=many_to_many_table_state,
            uid="ManyToMany1",
            record_uid=[
                [simple_table_state.id, ["simple_a_id"], {"uid": "SimpleRecord1"}],
                [simple_table_state.id, ["simple_b_id"], {"uid": "SimpleRecord2"}],
            ],
            record_data={
                FOREIGN_KEYS_LIST_COLUMN_NAME: [
                    [simple_table_state.id, ["simple_a_id"], {"uid": "SimpleRecord1"}],
                    [simple_table_state.id, ["simple_b_id"], {"uid": "SimpleRecord2"}],
                ],
            },
            event_type=WriteLogEventTypes.INSERT.value,
            created_at=datetime.datetime(
                2025, 5, 23, 12, 55, 00, tzinfo=datetime.timezone.utc
            ),
        ),
        WriteLog(
            migration_state=current_migration_state,
            table=many_to_many_table_state,
            uid="ManyToMany2",
            record_uid=[
                [simple_table_state.id, ["simple_a_id"], {"uid": "SimpleRecord2"}],
                [simple_table_state.id, ["simple_b_id"], {"uid": "SimpleRecord2"}],
            ],
            record_data={
                FOREIGN_KEYS_LIST_COLUMN_NAME: [
                    [simple_table_state.id, ["simple_a_id"], {"uid": "SimpleRecord2"}],
                    [simple_table_state.id, ["simple_b_id"], {"uid": "SimpleRecord2"}],
                ],
            },
            event_type=WriteLogEventTypes.INSERT.value,
            created_at=datetime.datetime(
                2025, 5, 23, 12, 56, 00, tzinfo=datetime.timezone.utc
            ),
        ),
        WriteLog(
            migration_state=current_migration_state,
            table=many_to_many_table_state,
            uid="ManyToMany3",
            record_uid=[
                [simple_table_state.id, ["simple_a_id"], {"uid": "SimpleRecord2"}],
                [simple_table_state.id, ["simple_b_id"], {"uid": "SimpleRecord2"}],
            ],
            record_data={
                FOREIGN_KEYS_LIST_COLUMN_NAME: [
                    [simple_table_state.id, ["simple_a_id"], {"uid": "SimpleRecord2"}],
                    [simple_table_state.id, ["simple_b_id"], {"uid": "SimpleRecord1"}],
                ],
            },
            event_type=WriteLogEventTypes.UPDATE.value,
            created_at=datetime.datetime(
                2025, 5, 23, 12, 56, 30, tzinfo=datetime.timezone.utc
            ),
        ),
        WriteLog(
            migration_state=current_migration_state,
            table=many_to_many_table_state,
            uid="ManyToMany4",
            record_uid=[
                [simple_table_state.id, ["simple_a_id"], {"uid": "SimpleRecord2"}],
                [simple_table_state.id, ["simple_b_id"], {"uid": "SimpleRecord1"}],
            ],
            record_data=None,
            event_type=WriteLogEventTypes.DELETE.value,
            created_at=datetime.datetime(
                2025, 5, 23, 12, 56, 45, tzinfo=datetime.timezone.utc
            ),
        ),
    ]

    for write_log_entry in input_write_log:
        replayer.replay(write_log_entry)
        write_log_entry.save()

    write_log = WriteLog.objects.all().order_by("id")

    assert len(write_log) == 6
    assert list(write_log) == input_write_log

    cursor.execute(f"SELECT id, uid FROM {simple_table}")  # noqa: S608

    rows = cursor.fetchall()

    assert [r[1] for r in rows] == ["SimpleRecord1", "SimpleRecord2"]

    simple_uid_to_id = {r[1]: r[0] for r in rows}

    cursor.execute(f"SELECT id, simple_a_id, simple_b_id FROM {many_to_many_table}")  # noqa: S608

    rows = cursor.fetchall()

    assert len(rows) == 1

    assert rows[0][1:] == (
        simple_uid_to_id["SimpleRecord1"],
        simple_uid_to_id["SimpleRecord2"],
    )
