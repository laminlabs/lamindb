from abc import ABC, abstractmethod
from copy import deepcopy

from django.apps import apps
from django.db.backends.utils import CursorWrapper
from django.db.models import ManyToManyField
from typing_extensions import override

from ._types import KeyConstraint, TableUID, UIDColumns


class DatabaseMetadataWrapper(ABC):
    """Provides an interface to metadata about the database's tables.

    Some information about the database can be retrieved from Django, but we
    also need to know things like key constraints and trigger information that
    are easier to retrieve from the database's information schema directly.
    This class provides a unified means for other history functions to access
    that metadata.
    """

    @abstractmethod
    def get_column_names(self, table: str, cursor: CursorWrapper) -> set[str]:
        raise NotImplementedError()

    @abstractmethod
    def get_tables_with_installed_triggers(self, cursor: CursorWrapper) -> set[str]:
        raise NotImplementedError()

    @abstractmethod
    def get_table_key_constraints(
        self, table: str, cursor: CursorWrapper
    ) -> tuple[KeyConstraint, list[KeyConstraint]]:
        raise NotImplementedError()

    def get_db_tables(self) -> set[str]:
        """Get the schema and table names for all tables in the database for which we might want to track history."""
        all_models = apps.get_models()

        tables: set[str] = set()

        for model in all_models:
            # Record the model's underlying table, as well as
            # the underlying table for all of its many-to-many
            # relationships
            tables.add(model._meta.db_table)

            for field in model._meta.many_to_many:
                m2m_field = field
                through_model = m2m_field.remote_field.through  # type: ignore

                tables.add(through_model._meta.db_table)

        return {t for t in tables if not t.startswith("django_")}

    def get_many_to_many_db_tables(self) -> set[str]:
        """Get the table names for all many-to-many tables in the database."""
        many_to_many_tables: set[str] = set()

        for model in apps.get_models():
            for field in model._meta.get_fields():
                if isinstance(field, ManyToManyField):
                    through_model = field.remote_field.through  # type: ignore
                    through_table = through_model._meta.db_table

                    many_to_many_tables.add(through_table)

        return many_to_many_tables

    def get_uid_columns(self, table: str, cursor: CursorWrapper) -> UIDColumns:
        """Get the UID columns for a given table."""
        if table in ("lamindb_paramvalue", "lamindb_featurevalue"):
            # Param and feature values should be uniquely identifiable by their value and hash
            return [
                TableUID(
                    source_table_name=table,
                    uid_columns=["value", "created_at"],
                    key_constraint=None,
                )
            ]
        else:
            column_names = self.get_column_names(table, cursor)

            # If the table has a 'uid' column, use that
            if "uid" in column_names:
                return [
                    TableUID(
                        source_table_name=table,
                        uid_columns=["uid"],
                        key_constraint=None,
                    )
                ]

            # Many-to-many tables are defined by the UIDs of the records pointed to by their
            # foreign-key constraints.

            many_to_many_tables = self.get_many_to_many_db_tables()

            uid_columns: UIDColumns = []

            if table in many_to_many_tables:
                _, foreign_key_constraints = self.get_table_key_constraints(
                    table=table, cursor=cursor
                )

                for constraint in foreign_key_constraints:
                    constraint_uid_columns = self.get_uid_columns(
                        table=constraint.target_table, cursor=cursor
                    )

                    if len(constraint_uid_columns) > 1:
                        raise ValueError(
                            "Many-to-many tables that reference other many-to-many tables aren't supported. "
                            f"'{table}' references '{constraint.target_table}', and both are many-to-many"
                        )

                    table_uid_for_constraint = deepcopy(constraint_uid_columns[0])
                    table_uid_for_constraint.key_constraint = constraint

                    uid_columns.append(table_uid_for_constraint)

                return uid_columns

        # If we hit this point, we don't know how to extract a unique identifier
        # for this table, so we need to panic.
        raise ValueError(f"Table '{table}' does not have known unique ID column(s)")


class PostgresDatabaseMetadataWrapper(DatabaseMetadataWrapper):
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
    CASE WHEN tc.contype = 'f' THEN af.attname ELSE NULL END AS target_column,
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
                source_column,
                target_column,
                target_table,
            ) = k

            if constraint_type == "PRIMARY KEY":
                if primary_key_constraint is None:
                    primary_key_constraint = KeyConstraint(
                        constraint_name=constraint_name,
                        constraint_type=constraint_type,
                        source_columns=[],
                        target_columns=[],
                        target_table=target_table,
                    )

                primary_key_constraint.source_columns.append(source_column)
                primary_key_constraint.target_columns.append(target_column)
            elif constraint_type == "FOREIGN KEY":
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

    @override
    def get_column_names(self, table: str, cursor: CursorWrapper) -> set[str]:
        cursor.execute(
            "SELECT column_name FROM information_schema.columns WHERE TABLE_NAME = %s ORDER BY ordinal_position",
            (table,),
        )

        return {r[0] for r in cursor.fetchall()}

    @override
    def get_tables_with_installed_triggers(self, cursor: CursorWrapper) -> set[str]:
        cursor.execute("""
SELECT
    DISTINCT event_object_table
FROM information_schema.triggers
WHERE trigger_name like 'lamindb_history_%';
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
    def get_column_names(self, table: str, cursor: CursorWrapper) -> set[str]:
        raise NotImplementedError()

    @override
    def get_tables_with_installed_triggers(self, cursor: CursorWrapper) -> set[str]:
        raise NotImplementedError()

    @override
    def get_table_key_constraints(
        self, table: str, cursor: CursorWrapper
    ) -> tuple[KeyConstraint, list[KeyConstraint]]:
        raise NotImplementedError()
