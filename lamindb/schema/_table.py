import lnschema_core

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
            if not hasattr(table, "__table__"):
                continue
            orm_to_model[table.__name__] = table
            tablename_to_orm[table.__table__.name] = table.__name__
            orm_to_tablename[table.__name__] = table.__table__.name

    @classmethod
    def convert_to_orm(cls, table_name: str):
        """Convert to camel case class name."""
        return (
            table_name
            if table_name.lower() != table_name
            else cls.tablename_to_orm.get(table_name)
        )

    @classmethod
    def get_model(cls, table_name: str):
        name = cls.convert_to_orm(table_name)
        if name is None:
            raise KeyError(f"Table {table_name} does NOT exist!")
        model = cls.orm_to_model.get(name)
        if model is None:
            raise AssertionError(f"Table {name} does NOT exist!")
        return model
