from abc import ABC, abstractmethod
from copy import deepcopy

from django.apps import apps
from django.db.backends.utils import CursorWrapper
from django.db.models import ManyToManyField
from typing_extensions import override

from ._types import Column, ColumnType, KeyConstraint, TableUID, UIDColumns


class DatabaseMetadataWrapper(ABC):
    """Provides an interface to metadata about the database's tables.

    Some information about the database can be retrieved from Django, but we
    also need to know things like key constraints and trigger information that
    are easier to retrieve from the database's information schema directly.
    This class provides a unified means for other write log functions to access
    that metadata.
    """

    def __init__(self):
        self._db_tables: set[str] | None = None
        self._many_to_many_tables: set[str] | None = None

    @abstractmethod
    def get_columns(self, table: str, cursor: CursorWrapper) -> set[Column]:
        raise NotImplementedError()

    def get_column_names(self, table: str, cursor: CursorWrapper) -> set[str]:
        return {c.name for c in self.get_columns(table, cursor)}

    @abstractmethod
    def get_tables_with_installed_triggers(self, cursor: CursorWrapper) -> set[str]:
        raise NotImplementedError()

    @abstractmethod
    def get_table_key_constraints(
        self, table: str, cursor: CursorWrapper
    ) -> tuple[KeyConstraint, list[KeyConstraint]]:
        raise NotImplementedError()

    def get_db_tables(self) -> set[str]:
        """Get the schema and table names for all tables in the database for which we might want to track writes."""
        if self._db_tables is None:
            all_models = apps.get_models()

            tables: set[str] = set()

            for model in all_models:
                # SQLRecord the model's underlying table, as well as
                # the underlying table for all of its many-to-many
                # relationships
                tables.add(model._meta.db_table)

                for field in model._meta.many_to_many:
                    m2m_field = field
                    through_model = m2m_field.remote_field.through  # type: ignore

                    tables.add(through_model._meta.db_table)

            self._db_tables = {t for t in tables if not t.startswith("django_")}

        return self._db_tables

    def get_many_to_many_db_tables(self) -> set[str]:
        """Get the table names for all many-to-many tables in the database."""
        if self._many_to_many_tables is None:
            many_to_many_tables: set[str] = set()

            for model in apps.get_models():
                for field in model._meta.get_fields():
                    if isinstance(field, ManyToManyField):
                        through_model = field.remote_field.through  # type: ignore
                        through_table = through_model._meta.db_table

                        many_to_many_tables.add(through_table)

            self._many_to_many_tables = many_to_many_tables

        return self._many_to_many_tables

    def _get_columns_by_name(
        self, table: str, column_names: list[str], cursor: CursorWrapper
    ) -> list[Column]:
        columns = self.get_columns(table=table, cursor=cursor)

        column_list: list[Column] = []

        for column_name in column_names:
            column = next((c for c in columns if c.name == column_name), None)

            if column is None:
                raise ValueError(
                    f"Table '{table}' doesn't have a column named '{column_name}'"
                )

            column_list.append(column)

        return column_list

    def get_uid_columns(self, table: str, cursor: CursorWrapper) -> UIDColumns:
        """Get the UID columns for a given table."""
        if table == "lamindb_featurevalue":
            _, foreign_key_constraints = self.get_table_key_constraints(
                table=table, cursor=cursor
            )

            for constraint in foreign_key_constraints:
                if constraint.target_table == "lamindb_feature":
                    feature_table_uid_columns = self.get_uid_columns(
                        table=constraint.target_table, cursor=cursor
                    )

                    if len(feature_table_uid_columns) > 1:
                        raise ValueError(
                            "Expected 'lamindb_feature' to have a simple UID, but got a complex one"
                        )

                    feature_table_uid = deepcopy(feature_table_uid_columns[0])
                    feature_table_uid.key_constraint = constraint

                    return [
                        feature_table_uid,
                        TableUID(
                            source_table_name=table,
                            uid_columns=self._get_columns_by_name(
                                table, ["hash"], cursor
                            ),
                            key_constraint=None,
                        ),
                    ]

            raise ValueError(
                "Expected lamindb_featurevalue to foreign-key to lamindb_feature, "
                "but found no such foreign key constraint"
            )
        else:
            column_names = self.get_column_names(table, cursor)

            # If the table has a 'uid' column, use that
            if "uid" in column_names:
                return [
                    TableUID(
                        source_table_name=table,
                        uid_columns=self._get_columns_by_name(table, ["uid"], cursor),
                        key_constraint=None,
                    )
                ]

            # Many-to-many tables are defined by the UIDs of the records pointed to by their
            # foreign-key constraints.

            many_to_many_tables = self.get_many_to_many_db_tables()
            many_to_many_tables.add("lamindb_recordjson")

            uid_columns: UIDColumns = []

            if table in many_to_many_tables:
                _, foreign_key_constraints = self.get_table_key_constraints(
                    table=table, cursor=cursor
                )

                for constraint in foreign_key_constraints:
                    if constraint.target_table in self.get_many_to_many_db_tables():
                        raise ValueError(
                            "Many-to-many tables that reference other many-to-many tables aren't supported. "
                            f"'{table}' references '{constraint.target_table}', and both are many-to-many"
                        )

                    constraint_uid_columns = self.get_uid_columns(
                        table=constraint.target_table, cursor=cursor
                    )

                    for uid_column in constraint_uid_columns:
                        table_uid_for_constraint = deepcopy(uid_column)
                        table_uid_for_constraint.key_constraint = constraint
                        uid_columns.append(table_uid_for_constraint)

                return uid_columns

        # If we hit this point, we don't know how to extract a unique identifier
        # for this table, so we need to panic.
        raise ValueError(f"Table '{table}' does not have known unique ID column(s)")


class PostgresDatabaseMetadataWrapper(DatabaseMetadataWrapper):
    def __init__(self) -> None:
        super().__init__()

        self._columns: dict[str, set[Column]] | None = None

    @override
    def get_table_key_constraints(
        self, table: str, cursor: CursorWrapper
    ) -> tuple[KeyConstraint, list[KeyConstraint]]:
        cursor.execute(
            """
SELECT
    tc.conname AS constraint_name,
    CASE tc.contype
        WHEN 'p' THEN 'PRIMARY KEY'
        WHEN 'f' THEN 'FOREIGN KEY'
    END AS constraint_type,
    a.attname AS source_column,
    a.attnum AS source_column_position,
    pg_catalog.format_type(a.atttypid, a.atttypmod) AS source_column_type,
    CASE WHEN tc.contype = 'f' THEN af.attname ELSE NULL END AS target_column,
    CASE WHEN tc.contype = 'f' THEN af.attnum ELSE NULL END AS target_column_position,
    CASE WHEN tc.contype = 'f' THEN pg_catalog.format_type(af.atttypid, af.atttypmod) ELSE NULL END AS target_column_type,
    CASE WHEN tc.contype = 'f' THEN tf.relname ELSE NULL END AS target_table
FROM
    pg_constraint tc
    JOIN pg_class t ON t.oid = tc.conrelid
    JOIN pg_namespace n ON n.oid = t.relnamespace
    JOIN pg_attribute a ON a.attrelid = tc.conrelid AND a.attnum = ANY(tc.conkey)
    LEFT JOIN pg_class tf ON tf.oid = tc.confrelid
    LEFT JOIN pg_attribute af ON af.attrelid = tc.confrelid AND
        af.attnum = tc.confkey[array_position(tc.conkey, a.attnum)]
WHERE
    t.relname = %s
    AND tc.contype IN ('p', 'f')
ORDER BY
    tc.conname,
    array_position(tc.conkey, a.attnum);
""",
            (table,),
        )

        keys = cursor.fetchall()

        primary_key_constraint = None

        foreign_key_constraints: dict[str, KeyConstraint] = {}

        for k in keys:
            (
                constraint_name,
                constraint_type,
                source_column_name,
                source_column_position,
                source_column_type,
                target_column_name,
                target_column_position,
                target_column_type,
                target_table,
            ) = k

            source_column = Column(
                name=source_column_name,
                type=self._get_column_type(source_column_type),
                ordinal_position=source_column_position,
            )

            if constraint_type == "PRIMARY KEY":
                if target_table is not None or target_column_name is not None:
                    raise Exception(
                        "Expected foreign key's target table/column to be NULL"
                    )

                if primary_key_constraint is None:
                    primary_key_constraint = KeyConstraint(
                        constraint_name=constraint_name,
                        constraint_type=constraint_type,
                        source_columns=[],
                        target_columns=[],
                        target_table=target_table,
                    )

                primary_key_constraint.source_columns.append(source_column)

            elif constraint_type == "FOREIGN KEY":
                target_column = Column(
                    name=target_column_name,
                    type=self._get_column_type(target_column_type),
                    ordinal_position=target_column_position,
                )

                if constraint_name not in foreign_key_constraints:
                    foreign_key_constraints[constraint_name] = KeyConstraint(
                        constraint_name=constraint_name,
                        constraint_type="FOREIGN KEY",
                        source_columns=[],
                        target_columns=[],
                        target_table=target_table,
                    )

                foreign_key_constraints[constraint_name].source_columns.append(
                    source_column
                )
                foreign_key_constraints[constraint_name].target_columns.append(
                    target_column
                )
            else:
                raise ValueError(f"Unhandled constraint type '{constraint_type}'")

        if primary_key_constraint is None:
            raise ValueError(
                f"Expected table {table} to have a primary key, but found none"
            )

        return (primary_key_constraint, list(foreign_key_constraints.values()))

    def _get_column_type(self, data_type: str) -> ColumnType:
        column_type: ColumnType

        if data_type in ("smallint", "integer", "bigint"):
            column_type = ColumnType.INT
        elif data_type in ("boolean",):
            column_type = ColumnType.BOOL
        elif data_type in ("character varying", "text"):
            column_type = ColumnType.STR
        elif data_type in ("jsonb",):
            column_type = ColumnType.JSON
        elif data_type in ("date",):
            column_type = ColumnType.DATE
        elif data_type in ("timestamp with time zone",):
            column_type = ColumnType.TIMESTAMPTZ
        elif data_type in ("double precision",):
            column_type = ColumnType.FLOAT
        else:
            raise ValueError(
                f"Don't know how to canonicalize column type '{data_type}'"
            )

        return column_type

    @override
    def get_columns(self, table: str, cursor: CursorWrapper) -> set[Column]:
        if self._columns is None:
            cursor.execute("""
SELECT
    table_name,
    column_name,
    data_type,
    ordinal_position
FROM information_schema.columns
WHERE
    table_schema not in ('pg_catalog', 'information_schema')
ORDER BY table_name, ordinal_position
""")
            self._columns = {}

            for row in cursor.fetchall():
                table_name, column_name, data_type, ordinal_position = row

                column_type = self._get_column_type(data_type)
                ordinal_position = ordinal_position

                if table_name not in self._columns:
                    self._columns[table_name] = set()

                self._columns[table_name].add(
                    Column(
                        name=column_name,
                        type=column_type,
                        ordinal_position=ordinal_position,
                    )
                )

        return self._columns[table]

    @override
    def get_tables_with_installed_triggers(self, cursor: CursorWrapper) -> set[str]:
        cursor.execute("""
SELECT
    DISTINCT event_object_table
FROM information_schema.triggers
WHERE trigger_name like 'lamindb_writelog_%';
""")
        tables = set()

        rows = cursor.fetchall()

        for row in rows:
            table_name = row[0]
            tables.add(table_name)

        return tables


class SQLiteDatabaseMetadataWrapper(DatabaseMetadataWrapper):
    """This is a placeholder until we implement the SQLite side of synchronization."""

    @override
    def get_columns(self, table: str, cursor: CursorWrapper) -> set[Column]:
        cursor.execute(
            """
SELECT
	name,
	"type",
	cid AS ordinal_position
FROM pragma_table_info(%s)
""",
            [table],
        )

        return {
            Column(
                name=r[0],
                type=self._data_type_to_column_type(r[1]),
                ordinal_position=r[2],
            )
            for r in cursor.fetchall()
        }

    @override
    def get_tables_with_installed_triggers(self, cursor: CursorWrapper) -> set[str]:
        cursor.execute("SELECT tbl_name FROM sqlite_master WHERE type = 'trigger';")

        return {r[0] for r in cursor.fetchall()}

    @override
    def get_table_key_constraints(
        self, table: str, cursor: CursorWrapper
    ) -> tuple[KeyConstraint, list[KeyConstraint]]:
        cursor.execute(
            """
SELECT
    p.name AS column_name,
    p.type AS data_type,
    p.cid AS ordinal_position
FROM sqlite_schema m
JOIN pragma_table_info(m.name) p ON p.pk > 0
WHERE m.type = 'table'
    AND m.name = %s
ORDER BY m.name, p.pk;
""",
            [table],
        )

        primary_key_columns = cursor.fetchall()

        primary_key_constraint = KeyConstraint(
            constraint_name="primary",
            constraint_type="PRIMARY KEY",
            source_columns=[
                Column(
                    name=row[0],
                    type=self._data_type_to_column_type(row[1]),
                    ordinal_position=row[2],
                )
                for row in primary_key_columns
            ],
            target_columns=[],
            target_table=table,
        )

        foreign_key_constraints: dict[int, KeyConstraint] = {}

        cursor.execute(
            """
SELECT
    id AS fk_id,
    "from" AS from_column,
    "table" AS referenced_table,
    "to" AS referenced_column
FROM pragma_foreign_key_list(%s)
ORDER BY id, seq;
""",
            [table],
        )

        rows = cursor.fetchall()

        for row in rows:
            fk_id, from_column, referenced_table, referenced_column = row

            if fk_id not in foreign_key_constraints:
                foreign_key_constraints[fk_id] = KeyConstraint(
                    constraint_name=f"foreign_key_{fk_id}",
                    constraint_type="FOREIGN KEY",
                    source_columns=[],
                    target_columns=[],
                    target_table=referenced_table,
                )

            foreign_key_constraints[fk_id].source_columns.append(from_column)
            foreign_key_constraints[fk_id].target_columns.append(referenced_column)

        return (
            primary_key_constraint,
            sorted(foreign_key_constraints.values(), key=lambda c: c.constraint_name),
        )

    def _data_type_to_column_type(self, data_type: str) -> ColumnType:
        if data_type.lower() in ("integer", "bigint", "smallint", "smallint unsigned"):
            return ColumnType.INT
        elif data_type.lower() in ("bool",):
            return ColumnType.BOOL
        elif data_type.lower().startswith("varchar") or data_type.lower() == "text":
            return ColumnType.STR
        elif data_type.lower() == "datetime":
            return ColumnType.TIMESTAMPTZ
        elif data_type.lower() == "date":
            return ColumnType.DATE
        elif data_type.lower() == "real":
            return ColumnType.FLOAT
        else:
            raise ValueError(f"Unhandled data type '{data_type}'")
