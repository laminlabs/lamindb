from types import ModuleType

from .. import schema


def _list_methods(module):
    return [getattr(module, i) for i in dir(module) if not i.startswith("__")]


alltables = {}
for schema_pkg in _list_methods(schema):
    alltables_pkg = _list_methods(schema_pkg)
    for table in alltables_pkg:
        if isinstance(table, (str, ModuleType)):
            continue
        alltables[table.__name__] = table
