from typing import Optional, Tuple

import lnschema_core
import sqlalchemy as sa
from lndb_setup import settings

from .. import schema


def _list_methods(module):
    return [getattr(module, i) for i in dir(module) if not i.startswith("_")]


class table_meta:
    orm_to_model = {}
    tablename_to_orm = {}
    orm_to_tablename = {}
    for schema_pkg in _list_methods(schema) + [lnschema_core]:
        alltables_pkg = _list_methods(schema_pkg)
        try:
            dev_tables = _list_methods(getattr(schema_pkg, "dev"))
        except AttributeError:
            dev_tables = []
        for table in alltables_pkg + dev_tables:
            if table.__class__.__name__ != "SQLModelMetaclass":
                continue
            orm_to_model[table.__name__] = table
            tablename_to_orm[table.__table__.name] = table.__name__
            orm_to_tablename[table.__name__] = table.__table__.name

    @classmethod
    def convert_to_tablename(cls, table_name: str):
        """Convert to lower case tablename."""
        return (
            table_name
            if table_name.lower() == table_name
            else cls.orm_to_tablename.get(table_name)
        )

    @classmethod
    def convert_to_orm(cls, table_name: str):
        """Convert to camel case class name."""
        return (
            table_name
            if table_name.lower() != table_name
            else cls.tablename_to_orm.get(table_name)
        )

    @classmethod
    def list_models(cls):
        return list(cls.orm_to_model.values())

    @classmethod
    def get_model(cls, table_name: str):
        name = cls.convert_to_orm(table_name)
        if name is None:
            raise KeyError(f"Table {table_name} does NOT exist!")
        model = cls.orm_to_model.get(name)
        if model is None:
            raise AssertionError(f"Table {name} does NOT exist!")
        return model

    @classmethod
    def get_pks(cls, table) -> list:
        if isinstance(table, str):
            model = cls.get_model(table_name=table)
        else:
            model = table

        return [i.name for i in model.__table__.primary_key.columns.values()]

    @classmethod
    def get_fields(cls, table):
        if isinstance(table, str):
            model = cls.get_model(table_name=table)
        else:
            model = table

        return list(model.__fields__.keys())

    @classmethod
    def get_foreign_keys(
        cls, table_name: str, inspector=None, referred: Tuple[str, str] = None
    ) -> dict:
        """Return foreign keys of a table.

        Returns {constrained_column: (referred_table, referred_column)}
        """
        table_name = cls.convert_to_tablename(table_name)

        if inspector is None:
            inspector = sa.inspect(settings.instance.db_engine())

        keys = {}
        results = inspector.get_foreign_keys(table_name)
        if len(results) > 0:
            for result in results:
                referred_table = result["referred_table"]
                for i, j in zip(
                    result["constrained_columns"], result["referred_columns"]
                ):
                    keys[i] = (referred_table, j)
        if referred is not None:
            referred_cols = (
                [referred[1]]
                if not isinstance(referred[1], (list, tuple))
                else list(referred[1])
            )
            keys = {
                result["constrained_columns"][0]: referred
                for result in results
                if (
                    result["referred_table"] == referred[0]
                    and result["referred_columns"] == referred_cols  # noqa
                )
            }
        return keys

    @classmethod
    def check_if_link_table(cls, table_name: str):
        """Check if a table is a link table.

        We define link tables there is overlap between primary and foreign keys
        """
        pks = cls.get_pks(table_name)
        fks = cls.get_foreign_keys(table_name).keys()
        intersect = set(pks).intersection(fks)
        if intersect:
            return intersect

    @classmethod
    def get_link_tables(cls, inspector=None) -> list:
        """Get all link tables."""
        if inspector is None:
            inspector = sa.inspect(settings.instance.db_engine())
        link_tables = []
        for name in inspector.get_table_names():
            if name == "nc_evolutions":  # table added by nocodb
                continue
            if cls.check_if_link_table(name):
                link_tables.append(name)

        return link_tables

    @classmethod
    def get_link_table(cls, table1: str, table2: str) -> Optional[str]:
        """Get the link table of two entities."""
        link_tables = cls.get_link_tables()
        pks = [f"{table1.split('.')[1]}_{i}" for i in cls.get_pks(table1)] + [
            f"{table2.split('.')[1]}_{i}" for i in cls.get_pks(table2)
        ]

        for table in link_tables:
            if set(pks) == set(cls.get_pks(table)):
                return table
        return None
