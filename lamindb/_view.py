import importlib
from typing import Optional

from IPython.display import display
from lamin_logger import colors
from lndb import settings
from lndb.dev._setup_schema import get_schema_module_name

from .dev.db._select import select


def view(n: int = 10, schema: Optional[str] = None):
    """View data.

    Args:
        n: Display the latest n rows of the tables.
        schema: Schema module to view. Default's to `None` and
            displays all schema modules.
    """
    if schema is not None:
        schema_names = [schema]
    else:
        schema_names = ["core"] + list(settings.instance.schema)

    for schema_name in schema_names:
        schema_module = importlib.import_module(get_schema_module_name(schema_name))

        tables = [
            i
            for i in schema_module.__dict__.values()
            if i.__class__.__name__ == "SQLModelMetaclass" and hasattr(i, "__table__")
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
