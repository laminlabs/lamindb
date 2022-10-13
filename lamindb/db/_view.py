from typing import Optional

from IPython.display import display
from lamin_logger import colors

from .. import schema
from ..dev.db._select import select
from ..schema import list_entities


def view(n: int = 10, schema_modules: Optional[list] = None):
    """View the latest edits to the database.

    Args:
        n: display the latest n rows of the tables
        schema_modules: a list of schema modules to view
            if None: will view all schema modules
    """
    if schema_modules is None:
        modules = [
            i
            for i in schema.__dir__()
            if getattr(schema, i).__class__.__name__ == "module"
            and not i.startswith("_")  # noqa
            and i != "core"  # noqa
        ]
        schema_modules = ["core"] + modules

    all_tables = list_entities()

    # core tables always display first
    for module in schema_modules:
        schema_module = getattr(schema, module)
        tables = [
            i
            for i in schema_module.__dir__()
            if getattr(schema_module, i).__class__.__name__ == "SQLModelMetaclass"
        ]
        if len(set(tables).intersection(all_tables)) == 0:
            continue
        section = f"* module: {colors.green(colors.bold(module))} *"
        section_no_color = f"* module: {module} *"
        print("*" * len(section_no_color))
        print(section)
        print("*" * len(section_no_color))
        for entity in tables:
            df = getattr(select, entity)().df()
            if df.shape[0] > 0:
                print(colors.blue(colors.bold(entity)))
                display(df.iloc[-n:])
