from .. import schema


def _list_methods(module):
    return [getattr(module, i) for i in dir(module) if not i.startswith("_")]


alltables = {}
for schema_pkg in _list_methods(schema):
    alltables_pkg = _list_methods(schema_pkg)
    for table in alltables_pkg:
        if table.__class__.__name__ != "SQLModelMetaclass":
            continue
        alltables[table.__name__] = table
