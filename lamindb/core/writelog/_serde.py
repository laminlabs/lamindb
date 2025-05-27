from datetime import datetime, timezone
from typing import Any

from lamindb.core.writelog._constants import FOREIGN_KEYS_LIST_COLUMN_NAME
from lamindb.core.writelog._db_metadata_wrapper import DatabaseMetadataWrapper
from lamindb.core.writelog._utils import get_latest_migration_state
from lamindb.models.writelog import WriteLog, WriteLogMigrationState, WriteLogTableState


class WriteLogDeserializationError(Exception):
    pass


class WriteLogSerDe:
    def __init__(
        self,
        db_metadata: DatabaseMetadataWrapper,
        table_id_mapping: dict[int, str] | None = None,
        latest_migration_state: WriteLogMigrationState | None = None,
    ):
        self._table_id_mapping = table_id_mapping
        self.db_metadata = db_metadata

        self._latest_migration_state = latest_migration_state

    @property
    def migration_state(self) -> WriteLogMigrationState:
        if self._latest_migration_state is None:
            self._latest_migration_state = get_latest_migration_state()

            if self._latest_migration_state is None:
                raise WriteLogDeserializationError(
                    "Write log migration state must be installed prior to deserializing write log entries"
                )

        return self._latest_migration_state

    @property
    def table_id_mapping(self) -> dict[int, str]:
        if self._table_id_mapping is None:
            self._table_id_mapping = {
                t.id: t.table_name for t in WriteLogTableState.objects.all()
            }

        return self._table_id_mapping

    def _remap_table(self, table_id: int) -> WriteLogTableState:
        try:
            table_name = self.table_id_mapping[table_id]
        except KeyError as e:
            raise KeyError(
                f"Unable to remap table ID {table_id} from "
                "the source, since that ID is not in the source's table state"
            ) from e

        try:
            table = WriteLogTableState.objects.get(table_name=table_name)
        except WriteLogTableState.DoesNotExist as e:
            raise WriteLogDeserializationError(
                f"Unable to locate table '{table_name}', which is "
                "present in the source, in this instance. Perhaps your "
                "schema has diverged with that of the source?"
            ) from e

        return table

    def _remap_foreign_key_uid(self, foreign_key_uid):
        table_id, column_names, uid = foreign_key_uid
        remapped_table_id = self._remap_table(table_id).id

        return [remapped_table_id, column_names, uid]

    def serialize(self, write_log_entry: WriteLog) -> dict[str, Any]:
        return {
            "migration_state_id": write_log_entry.migration_state.migration_state_id,
            "table_id": write_log_entry.table.id,
            "uid": write_log_entry.uid,
            "space_id": write_log_entry.space_id,
            "created_by_uid": write_log_entry.created_by_uid,
            "branch_id": write_log_entry.branch_id,
            "run_uid": write_log_entry.run_uid,
            "record_uid": write_log_entry.record_uid,
            "record_data": write_log_entry.record_data,
            "event_type": write_log_entry.event_type,
            "created_at": write_log_entry.created_at.astimezone(
                timezone.utc
            ).isoformat(),
        }

    def deserialize(self, serialized_write_log_entry: dict[str, Any]) -> WriteLog:
        try:
            entry_uid = serialized_write_log_entry["uid"]
        except KeyError as e:
            raise WriteLogDeserializationError(
                f"Serialized write log entry has no UID: '{serialized_write_log_entry}'"
            ) from e

        try:
            if (
                self.migration_state.migration_state_id
                != serialized_write_log_entry["migration_state_id"]
            ):
                raise WriteLogDeserializationError(
                    f"Write log entry '{entry_uid}' has a different migration state than this instance's "
                    "current migration state. This indicates that the schemas for the source and destination "
                    "have diverged. Migrating write logs across divergent schemas is not currently supported."
                )

            table = self._remap_table(serialized_write_log_entry["table_id"])

            if table.table_name in self.db_metadata.get_many_to_many_db_tables():
                record_uid = []

                for many_to_many_uid in serialized_write_log_entry["record_uid"]:
                    record_uid.append(self._remap_foreign_key_uid(many_to_many_uid))
            else:
                record_uid = serialized_write_log_entry["record_uid"]

            record_data = serialized_write_log_entry["record_data"]
            record_data[FOREIGN_KEYS_LIST_COLUMN_NAME] = [
                self._remap_foreign_key_uid(k)
                for k in record_data[FOREIGN_KEYS_LIST_COLUMN_NAME]
            ]

            return WriteLog(
                table=table,
                migration_state=self.migration_state,
                uid=entry_uid,
                space_id=serialized_write_log_entry["space_id"],
                created_by_uid=serialized_write_log_entry["created_by_uid"],
                branch_id=serialized_write_log_entry["branch_id"],
                run_uid=serialized_write_log_entry["run_uid"],
                record_uid=record_uid,
                record_data=record_data,
                event_type=serialized_write_log_entry["event_type"],
                created_at=datetime.fromisoformat(
                    serialized_write_log_entry["created_at"]
                ),
            )

        except KeyError as e:
            raise WriteLogDeserializationError(
                f"Serialized write log entry '{entry_uid}' is missing key '{str(e)}'"
            ) from e
