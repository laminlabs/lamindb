import importlib
import inspect
from typing import Optional

from IPython.display import display
from lamin_logger import colors
from lamindb_setup import settings
from lamindb_setup.dev._setup_schema import get_schema_module_name
from lnschema_core import BaseORM

from ._select import select


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

        orms = [
            i
            for i in schema_module.__dict__.values()
            if inspect.isclass(i) and issubclass(i, BaseORM) and i.__name__ != "BaseORM"
        ]
        if len(schema_names) > 1:
            section = f"* module: {colors.green(colors.bold(schema_name))} *"
            section_no_color = f"* module: {schema_name} *"
            print("*" * len(section_no_color))
            print(section)
            print("*" * len(section_no_color))
        for orm in orms:
            df = select(orm).df()
            if df.shape[0] > 0:
                print(colors.blue(colors.bold(orm.__name__)))
                display(df.iloc[-n:])
