from datetime import datetime, timezone
from typing import Generator

import pytest
from lamindb.core.writelog._constants import FOREIGN_KEYS_LIST_COLUMN_NAME
from lamindb.core.writelog._serde import WriteLogSerDe
from lamindb.core.writelog._trigger_installer import WriteLogEventTypes
from lamindb.models.writelog import (
    DEFAULT_BRANCH_CODE,
    DEFAULT_CREATED_BY_UID,
    DEFAULT_RUN_UID,
    WriteLog,
    WriteLogMigrationState,
    WriteLogTableState,
)

from .writelog_test_utils import FakeMetadataWrapper


@pytest.fixture(scope="function")
def migration_state() -> Generator[WriteLogMigrationState, None, None]:
    migration_state = WriteLogMigrationState.objects.create(
        migration_state_id={"test_app": 42}
    )

    yield migration_state

    migration_state.delete()


@pytest.fixture(scope="function")
def simple_table() -> Generator[WriteLogTableState, None, None]:
    table = WriteLogTableState.objects.create(
        table_name="simple_table", backfilled=True
    )

    yield table

    table.delete()


@pytest.fixture(scope="function")
def many_to_many_table() -> Generator[WriteLogTableState, None, None]:
    table = WriteLogTableState.objects.create(
        table_name="many_to_many_table", backfilled=True
    )

    yield table

    table.delete()


def assert_write_logs_equal(a: WriteLog, b: WriteLog):
    assert (
        a.migration_state.id == b.migration_state.id
        and a.table.id == b.table.id
        and a.uid == b.uid
        and a.space_uid == b.space_uid
        and a.created_by_uid == b.created_by_uid
        and a.branch_code == b.branch_code
        and a.run_uid == b.run_uid
        and a.record_uid == b.record_uid
        and a.record_data == b.record_data
        and a.event_type == b.event_type
        and a.created_at == b.created_at
    )


def test_local_writelog_simple_table(simple_table, migration_state):
    db_metadata = FakeMetadataWrapper()
    db_metadata.set_db_tables({simple_table.table_name})

    serde = WriteLogSerDe(
        db_metadata=db_metadata,
        table_id_mapping={simple_table.id: simple_table.table_name},
        latest_migration_state=migration_state,
    )

    write_log_entry = WriteLog(
        migration_state=migration_state,
        table=simple_table,
        uid="WriteLog1",
        record_uid=["writeloguid"],
        record_data={
            "int_col": 42,
            "str_col": "hello world",
            FOREIGN_KEYS_LIST_COLUMN_NAME: [],
        },
        event_type=WriteLogEventTypes.INSERT.value,
        created_at=datetime(2025, 5, 26, 12, 34, 56, tzinfo=timezone.utc),
    )

    serialized_entry = serde.serialize(write_log_entry)

    deserialized_entry = serde.deserialize(serialized_entry)

    assert_write_logs_equal(write_log_entry, deserialized_entry)


def test_local_writelog_many_to_many_table(
    many_to_many_table, simple_table, migration_state
):
    db_metadata = FakeMetadataWrapper()

    db_metadata.set_db_tables({simple_table.table_name, many_to_many_table.table_name})
    db_metadata.set_many_to_many_db_tables({many_to_many_table.table_name})

    serde = WriteLogSerDe(
        db_metadata=db_metadata,
        table_id_mapping={
            t.id: t.table_name for t in (simple_table, many_to_many_table)
        },
        latest_migration_state=migration_state,
    )

    write_log_entry = WriteLog(
        migration_state=migration_state,
        table=many_to_many_table,
        uid="WriteLog1",
        record_uid=[
            [simple_table.id, ["sample_ref_a"], {"uid": "foo"}],
            [simple_table.id, ["sample_ref_b"], {"uid": "bar"}],
        ],
        record_data={
            FOREIGN_KEYS_LIST_COLUMN_NAME: [
                [simple_table.id, ["sample_ref_a"], {"uid": "foo"}],
                [simple_table.id, ["sample_ref_b"], {"uid": "bar"}],
            ]
        },
        event_type=WriteLogEventTypes.INSERT.value,
        created_at=datetime(2025, 5, 26, 12, 34, 56, tzinfo=timezone.utc),
    )

    serialized_entry = serde.serialize(write_log_entry)

    deserialized_entry = serde.deserialize(serialized_entry)

    assert_write_logs_equal(write_log_entry, deserialized_entry)


def test_table_remapping(simple_table, many_to_many_table, migration_state):
    db_metadata = FakeMetadataWrapper()

    db_metadata.set_db_tables({simple_table.table_name, many_to_many_table.table_name})
    db_metadata.set_many_to_many_db_tables({many_to_many_table.table_name})

    REMOTE_SIMPLE_TABLE_ID = 8675309
    REMOTE_MANY_TO_MANY_TABLE_ID = 5008000

    serde = WriteLogSerDe(
        db_metadata=db_metadata,
        table_id_mapping={
            REMOTE_SIMPLE_TABLE_ID: simple_table.table_name,
            REMOTE_MANY_TO_MANY_TABLE_ID: many_to_many_table.table_name,
        },
        latest_migration_state=migration_state,
    )

    serialized_entry = {
        "migration_state_id": migration_state.migration_state_id,
        "table_id": REMOTE_MANY_TO_MANY_TABLE_ID,
        "uid": "RemoteLog1",
        "record_uid": [
            [REMOTE_SIMPLE_TABLE_ID, ["sample_ref_a"], {"uid": "foo"}],
            [REMOTE_SIMPLE_TABLE_ID, ["sample_ref_b"], {"uid": "bar"}],
        ],
        "record_data": {
            FOREIGN_KEYS_LIST_COLUMN_NAME: [
                [REMOTE_SIMPLE_TABLE_ID, ["sample_ref_a"], {"uid": "foo"}],
                [REMOTE_SIMPLE_TABLE_ID, ["sample_ref_b"], {"uid": "bar"}],
            ]
        },
        "event_type": WriteLogEventTypes.INSERT.value,
        "space_uid": None,
        "created_by_uid": DEFAULT_CREATED_BY_UID,
        "branch_code": DEFAULT_BRANCH_CODE,
        "run_uid": DEFAULT_RUN_UID,
        "created_at": "2025-05-26T12:34:56.000000+00:00",
    }

    expected_write_log_entry = WriteLog(
        migration_state=migration_state,
        table=many_to_many_table,
        uid="RemoteLog1",
        record_uid=[
            [simple_table.id, ["sample_ref_a"], {"uid": "foo"}],
            [simple_table.id, ["sample_ref_b"], {"uid": "bar"}],
        ],
        record_data={
            FOREIGN_KEYS_LIST_COLUMN_NAME: [
                [simple_table.id, ["sample_ref_a"], {"uid": "foo"}],
                [simple_table.id, ["sample_ref_b"], {"uid": "bar"}],
            ],
        },
        event_type=WriteLogEventTypes.INSERT.value,
        created_at=datetime(2025, 5, 26, 12, 34, 56, tzinfo=timezone.utc),
    )

    deserialized_entry = serde.deserialize(serialized_entry)

    assert_write_logs_equal(expected_write_log_entry, deserialized_entry)
