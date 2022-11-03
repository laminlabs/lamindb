import importlib
from typing import Optional

from IPython.display import display
from lamin_logger import colors
from lndb_setup import settings
from lndb_setup._setup_schema import get_schema_module_name

from .dev.db._select import select


def view(n: int = 10, schema_names: Optional[list] = None):
    """View the latest edits to the database.

    Args:
        n: display the latest n rows of the tables
        schema_names: a list of schema modules to view
            if None: will view all schema modules
    """
    if schema_names is None:
        if settings.instance.schema_modules is not None:
            schema_names = settings.instance.schema_modules.split(", ")
        else:
            schema_names = []

    for schema_name in ["core"] + schema_names:
        schema_module = importlib.import_module(get_schema_module_name(schema_name))

        tables = [
            i
            for i in schema_module.__dict__.values()
            if i.__class__.__name__ == "SQLModelMetaclass"
        ]
        section = f"* module: {colors.green(colors.bold(schema_name))} *"
        section_no_color = f"* module: {schema_name} *"
        print("*" * len(section_no_color))
        print(section)
        print("*" * len(section_no_color))
        for entity in tables:
            df = select(entity).df()
            if df.shape[0] > 0:
                print(colors.blue(colors.bold(entity.__name__)))
                display(df.iloc[-n:])
