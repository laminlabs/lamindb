from dataclasses import dataclass
from typing import Any

from django.db.backends.utils import CursorWrapper

from lamindb.core.writelog._constants import FOREIGN_KEYS_LIST_COLUMN_NAME
from lamindb.core.writelog._trigger_installer import WriteLogEventTypes
from lamindb.models.writelog import WriteLog, WriteLogTableState

from ._db_metadata_wrapper import DatabaseMetadataWrapper


@dataclass
class WriteLogForeignKey:
    table_name: str
    foreign_key_columns: list[str]
    foreign_uid: dict[str, str]


class WriteLogReplayer:
    def __init__(
        self,
        db_metadata: DatabaseMetadataWrapper,
        table_id_mapping: dict[int, int],
        cursor: CursorWrapper,
    ):
        self.db_metadata = db_metadata
        self.table_id_mapping = table_id_mapping
        self.cursor = cursor

    def replay(self, write_log_record: WriteLog):
        table_name = write_log_record.table.table_name

        resolved_record_uid = self._resolve_uid(table_name, write_log_record.record_uid)

        # If we're deleting the record, its (resolved) UID is the only thing we need
        if write_log_record.event_type == WriteLogEventTypes.DELETE:
            self.cursor.execute(f"""
                DELETE FROM {table_name}
                WHERE {" AND ".join("k = v" for k, v in resolved_record_uid.items())}
                LIMIT 1
            """)  # noqa: S608
        else:
            if write_log_record.record_data is None:
                raise ValueError(
                    f"Expected non-null record data for write log record {write_log_record.uid} (type: {write_log_record.event_type})"
                )

            record_data = self._build_record_data(write_log_record.record_data)

            # FIXME need DB-specific casts here too
            if write_log_record.event_type == WriteLogEventTypes.INSERT:
                record_keys = [d[0] for d in record_data]
                record_values = [d[1] for d in record_data]

                self.cursor.execute(f"""
                INSERT INTO {table_name}
                ({", ".join(record_keys)}) VALUES ({", ".join(record_values)})
                """)  # noqa: S608
            elif write_log_record.event_type == WriteLogEventTypes.UPDATE:
                set_statements = ", ".join(f"{k} = {v}" for (k, v) in record_data)

                self.cursor.execute(f"""
                UPDATE {table_name}
                SET {set_statements}
                WHERE {" AND ".join("k = v" for k, v in resolved_record_uid.items())} LIMIT 1
                """)  # noqa: S608
            else:
                raise ValueError(
                    f"Unhandled record event type {write_log_record.event_type}"
                )

    def _resolve_uid(self, table_name: str, record_uid) -> dict[str, int]:
        if table_name in self.db_metadata.get_many_to_many_db_tables():
            resolved_record_uid = self._resolve_foreign_keys_list(record_uid)
        else:
            table_uid_list = self.db_metadata.get_uid_columns(
                table=table_name, cursor=self.cursor
            )

            if (
                len(table_uid_list) != 1
                and table_uid_list[0].source_table_name == table_name
            ):
                raise ValueError(
                    f"Expected standard table {table_name}'s UID to refer only to itself"
                )

            uid_columns = table_uid_list[0].uid_columns

            if not (
                isinstance(record_uid, list) and len(record_uid) == len(uid_columns)
            ):
                raise ValueError(
                    f"Expected standard table {table_name}'s UID to be a list of length {len(uid_columns)}"
                )

            # FIXME need a DB-specific cast here
            lookup_where_clause = " AND ".join(
                f"k = {v}" for (k, v) in zip(uid_columns, record_uid)
            )

            primary_key, _ = self.db_metadata.get_table_key_constraints(
                table=table_name, cursor=self.cursor
            )

            self.cursor.execute(f"""
            SELECT {", ".join(primary_key.source_columns)}
            FROM {table_name}
            WHERE {lookup_where_clause}
            """)  # noqa: S608

        return resolved_record_uid

    def _resolve_foreign_keys_list(self, foreign_keys_list) -> dict[str, int]:
        """Resolves a reference to a foreign record by UID into the foreign-keys needed by the source table.

        Write logs refer to records in other tables by their UID. These references are stored as a list of 3-tuples,
        where each 3-tuple is:

            [
                ID of the destination table in WriteLogTableState,
                list of table fields containing the foreign key,
                record/value pairs defining the foreign record's UID
            ]

        This method resolves these three-tuples into a dict that maps column names to record IDs.

        This assumes that all tables are keyed with integers for simplicity.
        """
        if not isinstance(foreign_keys_list, list) or any(
            not (isinstance(x, list) and len(x) == 3) for x in foreign_keys_list
        ):
            raise ValueError("Expected a foreign keys list to be a list of 3-tuples")

        resolved_foreign_keys: dict[str, int] = {}

        for foreign_table_uid in foreign_keys_list:
            resolved_foreign_key = self._foreign_key_from_json(
                json_obj=foreign_table_uid
            )
            resolved_foreign_keys.update(
                self._get_foreign_key_values(resolved_foreign_key)
            )

        return resolved_foreign_keys

    def _build_record_data(self, record_data: dict[str, Any]) -> list[tuple[str, str]]:
        # We're outputting this data as tuples so that it's easy to output keys and values
        # as two separate lists with the same order in INSERT.
        record_data_tuples = []

        for key, value in record_data.items():
            if key == FOREIGN_KEYS_LIST_COLUMN_NAME:
                for (
                    foreign_key_column,
                    foreign_key_value,
                ) in self._resolve_foreign_keys_list(value).items():
                    record_data_tuples.append(
                        (foreign_key_column, str(foreign_key_value))
                    )
            else:
                record_data_tuples.append((key, value))

        return record_data_tuples

    def _foreign_key_from_json(self, json_obj: list) -> "WriteLogForeignKey":
        try:
            foreign_table_id_at_source, foreign_key_columns, foreign_uid = json_obj
            table_id = self.table_id_mapping[foreign_table_id_at_source]

            table_name = WriteLogTableState.objects.get(id=table_id).table_name

            return WriteLogForeignKey(
                table_name=table_name,
                foreign_key_columns=foreign_key_columns,
                foreign_uid=foreign_uid,
            )
        except ValueError as e:
            raise ValueError(
                f"Malformed write log foreign key '{json_obj}' encountered"
            ) from e

    def _get_foreign_key_values(
        self, foreign_key: WriteLogForeignKey
    ) -> dict[str, int]:
        primary_key, _ = self.db_metadata.get_table_key_constraints(
            foreign_key.table_name, cursor=self.cursor
        )

        uid_constraints = {f"{k} = '{v}'" for (k, v) in foreign_key.foreign_uid.items()}

        self.cursor.execute(
            f"SELECT {','.join(primary_key.source_columns)} "  # noqa: S608
            f"FROM {foreign_key.table_name} "
            f"WHERE {uid_constraints}"
        )

        rows = self.cursor.fetchall()

        if len(rows) == 0:
            raise ValueError(
                f"Record matching foreign key constraint {foreign_key} not found"
            )

        if len(rows) > 1:
            raise ValueError(
                f"Found more than one record matching foreign key constraint {foreign_key}"
            )

        row = rows[0]
        return dict(zip(primary_key.source_columns, row))

    def replay_standard(self, write_log_record: WriteLog):
        pass
