import datetime
from dataclasses import dataclass
from typing import Any

from django.db import connection
from django.db.backends.utils import CursorWrapper

from lamindb.core.writelog._constants import FOREIGN_KEYS_LIST_COLUMN_NAME
from lamindb.core.writelog._trigger_installer import WriteLogEventTypes
from lamindb.core.writelog._types import Column, ColumnType
from lamindb.models.writelog import TableState, WriteLog

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
        cursor: CursorWrapper,
    ):
        self.db_metadata = db_metadata
        self.cursor = cursor

    def replay(self, write_log_record: WriteLog):
        table_name = write_log_record.table.table_name

        # If we're deleting the record, its (resolved) UID is the only thing we need
        if write_log_record.event_type == WriteLogEventTypes.DELETE.value:
            record_id = self._resolve_uid(table_name, write_log_record.record_uid)

            self.cursor.execute(f"""
                DELETE FROM {table_name}
                WHERE {" AND ".join(f"{k} = {v}" for k, v in record_id.items())}
                LIMIT 1
            """)  # noqa: S608
        else:
            record_data = self._build_record_data(write_log_record)

            if write_log_record.event_type == WriteLogEventTypes.INSERT.value:
                record_column_names = [d[0].name for d in record_data]
                record_values = [
                    self._cast_column(column=d[0], value=d[1]) for d in record_data
                ]

                if table_name in self.db_metadata.get_many_to_many_db_tables():
                    record_id = self._resolve_foreign_keys_list(
                        write_log_record.record_uid
                    )

                    for column_name, value in record_id.items():
                        record_column_names.append(column_name)
                        record_values.append(str(value))
                else:
                    uid_list = self.db_metadata.get_uid_columns(
                        table=table_name, cursor=self.cursor
                    )

                    if len(uid_list) != 1:
                        raise ValueError(
                            f"Table {table_name} is not marked as many-to-many, but has "
                            "more than one table reference in its UID"
                        )

                    record_uid = uid_list[0]

                    if len(record_uid.uid_columns) != len(write_log_record.record_uid):
                        raise ValueError(
                            f"Write log record {write_log_record.uid} expected to "
                            f"have {len(record_uid.uid_columns)} components to its "
                            "UID, but only had {len(write_log_record.record_uid)}"
                        )

                    for uid_column, value in zip(
                        record_uid.uid_columns, write_log_record.record_uid
                    ):
                        record_column_names.append(uid_column.name)
                        record_values.append(self._cast_column(uid_column, value))

                self.cursor.execute(f"""
                INSERT INTO {table_name}
                ({", ".join(record_column_names)}) VALUES ({", ".join(record_values)})
                """)  # noqa: S608
            elif write_log_record.event_type == WriteLogEventTypes.UPDATE.value:
                record_id = self._resolve_uid(table_name, write_log_record.record_uid)

                set_statements = ", ".join(
                    f"{column.name} = {self._cast_column(column, value)}"
                    for (column, value) in record_data
                )

                self.cursor.execute(f"""
                UPDATE {table_name}
                SET {set_statements}
                WHERE {" AND ".join(f"{k} = {v}" for (k, v) in record_id.items())} LIMIT 1
                """)  # noqa: S608
            else:
                raise ValueError(
                    f"Unhandled record event type {write_log_record.event_type}"
                )

    def _cast_column(self, column: Column, value: Any) -> str:
        """Returns a string representation of the column that will cast it to the appropriate type."""
        if connection.vendor == "postgresql":
            return self._cast_column_postgres(column, value)
        elif connection.vendor == "sqlite":
            return self._cast_column_sqlite(column, value)
        else:
            raise ValueError(f"Unsupported connection vendor '{connection.vendor}'")

    def _cast_column_postgres(self, column: Column, value: Any) -> str:
        if value is None:
            return "NULL"
        elif column.type == ColumnType.INT:
            return str(value)
        elif column.type == ColumnType.BOOL:
            return "TRUE" if value is True else "FALSE"
        elif column.type == ColumnType.STR:
            return f"'{value}'"
        elif column.type == ColumnType.DATE:
            return f"date('{value}')"
        elif column.type == ColumnType.FLOAT:
            return str(value)
        elif column.type == ColumnType.JSON:
            return f"to_jsonb('{value}'::text)"
        elif column.type == ColumnType.TIMESTAMPTZ:
            return f"timestamptz('{value}')"
        else:
            raise ValueError(f"Unhandled type {column.type}")

    def _cast_column_sqlite(self, column: Column, value: Any) -> str:
        if value is None:
            return "NULL"
        elif column.type == ColumnType.INT:
            return str(value)
        elif column.type == ColumnType.BOOL:
            return "1" if value is True else "0"
        elif column.type == ColumnType.STR:
            return f"'{value}'"
        elif column.type == ColumnType.DATE:
            return f"'{value}'"
        elif column.type == ColumnType.FLOAT:
            return f"CAST('{str(value)}' AS REAL)"
        elif column.type == ColumnType.JSON:
            return f"'{value}'"
        elif column.type == ColumnType.TIMESTAMPTZ:
            formatted_datetime = (
                datetime.datetime.fromisoformat(value)
                .astimezone(datetime.timezone.utc)
                .strftime("%Y-%m-%d %H:%M:%S.%f")
            )

            return f"'{formatted_datetime}'"
        else:
            raise ValueError(f"Unhandled type {column.type}")

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

            lookup_where_clause = " AND ".join(
                f"{k.name} = {self._cast_column(k, v)}"
                for (k, v) in zip(uid_columns, record_uid)
            )

            primary_key, _ = self.db_metadata.get_table_key_constraints(
                table=table_name, cursor=self.cursor
            )

            self.cursor.execute(f"""
            SELECT {", ".join(c.name for c in primary_key.source_columns)}
            FROM {table_name}
            WHERE {lookup_where_clause}
            LIMIT 1
            """)  # noqa: S608

            table_primary_key_values = self.cursor.fetchone()

            if table_primary_key_values is None:
                raise ValueError(
                    f"Unable to locate a record in {table_name} with UID {record_uid}"
                )

            return dict(
                zip(
                    [c.name for c in primary_key.source_columns],
                    table_primary_key_values,
                )
            )

        return resolved_record_uid

    def _build_record_data(self, record: WriteLog) -> list[tuple[Column, str]]:
        record_data: dict[str, Any] | None = record.record_data

        if record_data is None:
            raise ValueError(
                f"Expected non-null record data for write log record {record.uid} (type: {record.event_type})"
            )

        # We're outputting this data as tuples so that it's easy to output keys and values
        # as two separate lists with the same order in INSERT.
        record_data_tuples: list[tuple[Column, str]] = []

        columns_by_name = {
            c.name: c
            for c in self.db_metadata.get_columns(
                table=record.table.table_name, cursor=self.cursor
            )
        }

        for key, value in record_data.items():
            if key == FOREIGN_KEYS_LIST_COLUMN_NAME:
                for (
                    foreign_key_column,
                    foreign_key_value,
                ) in self._resolve_foreign_keys_list(value).items():
                    if foreign_key_column not in columns_by_name:
                        raise ValueError(
                            f"Table {record.table.table_name} does not have a column named {foreign_key_column}"
                        )

                    record_data_tuples.append(
                        (columns_by_name[foreign_key_column], str(foreign_key_value))
                    )
            else:
                if key not in columns_by_name:
                    raise ValueError(
                        f"Table {record.table.table_name} does not have a column named {key}"
                    )

                record_data_tuples.append((columns_by_name[key], value))

        return record_data_tuples

    def _resolve_foreign_keys_list(self, foreign_keys_list) -> dict[str, int]:
        """Resolves a reference to a foreign record by UID into the foreign-keys needed by the source table.

        Write logs refer to records in other tables by their UID. These references are stored as a list of 3-tuples,
        where each 3-tuple is:

            [
                ID of the destination table in TableState,
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

    def _foreign_key_from_json(self, json_obj: list) -> "WriteLogForeignKey":
        try:
            table_id, foreign_key_columns, foreign_uid = json_obj

            if not isinstance(table_id, int):
                raise ValueError(
                    f"Expected the first element of a foreign key's JSON representation to be a table state ID (JSON: {json_obj})"
                )

            if not isinstance(foreign_key_columns, list):
                raise ValueError(
                    f"Expected the second element of a foreign key's JSON representation to be a list of column names (JSON: {json_obj})"
                )

            if not isinstance(foreign_uid, dict):
                raise ValueError(
                    f"Expected the third element of a foreign key's JSON representation to be a dict mapping columns to values (JSON: {json_obj})"
                )

            table_name = TableState.objects.get(id=table_id).table_name

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

        if len(foreign_key.foreign_key_columns) != len(primary_key.source_columns):
            raise ValueError(
                f"Expected number of primary key columns for {primary_key.target_table} "
                f"({len(primary_key.source_columns)}) and number of foreign key columns "
                f"for {foreign_key.table_name} ({len(foreign_key.foreign_key_columns)}) to match"
            )

        uid_constraints = " AND ".join(
            f"{k} = '{v}'" for (k, v) in foreign_key.foreign_uid.items()
        )

        self.cursor.execute(
            f"SELECT {','.join(c.name for c in primary_key.source_columns)} "  # noqa: S608
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

        return dict(zip(foreign_key.foreign_key_columns, row))
