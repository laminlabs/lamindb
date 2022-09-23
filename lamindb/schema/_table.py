from .. import schema


def _list_methods(module):
    return [getattr(module, i) for i in dir(module) if not i.startswith("_")]


class Table:
    all = {}
    for schema_pkg in _list_methods(schema):
        alltables_pkg = _list_methods(schema_pkg)
        for table in alltables_pkg:
            if table.__class__.__name__ != "SQLModelMetaclass":
                continue
            all[table.__name__] = table

    @classmethod
    def list_models(cls):
        return list(cls.all.values())

    @classmethod
    def get_model(cls, table_name):
        model = cls.all.get(table_name)
        if model is None:
            raise AssertionError(f"Table {table_name} does NOT exist!")
        return model

    @classmethod
    def get_pks(cls, table):
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
