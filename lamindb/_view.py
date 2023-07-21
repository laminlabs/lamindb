import importlib
import inspect
from typing import List, Optional

from IPython.display import display
from lamin_utils import colors
from lamindb_setup import settings
from lamindb_setup.dev._setup_schema import get_schema_module_name
from lnschema_core import ORM

from ._select import select


def view(n: int = 10, schema: Optional[str] = None, orms: Optional[List[str]] = None):
    """View data.

    Args:
        n: ``int = 10`` Display the last `n` rows of a table.
        schema: ``Optional[str] = None`` Schema module to view. Default's to
            `None` and displays all schema modules.
        orms: ``Optional[List[str]] = None`` List of ORM names. Defaults to
            `None` and lists all ORMs.

    Examples:
        >>> ln.view()
    """
    if schema is not None:
        schema_names = [schema]
    else:
        schema_names = ["core"] + list(settings.instance.schema)

    for schema_name in schema_names:
        schema_module = importlib.import_module(get_schema_module_name(schema_name))

        all_orms = {
            orm
            for orm in schema_module.__dict__.values()
            if inspect.isclass(orm) and issubclass(orm, ORM) and orm.__name__ != "ORM"
        }
        if orms is not None:
            filtered_orms = {orm for orm in all_orms if orm.__name__ in orms}
        else:
            filtered_orms = all_orms
        if len(schema_names) > 1:
            section = f"* module: {colors.green(colors.bold(schema_name))} *"
            section_no_color = f"* module: {schema_name} *"
            print("*" * len(section_no_color))
            print(section)
            print("*" * len(section_no_color))
        for orm in sorted(filtered_orms, key=lambda x: x.__name__):
            if hasattr(orm, "updated_at"):
                df = select(orm).order_by("-updated_at")[:n].df()
            else:
                # need to adjust in the future
                df = select(orm).df().iloc[-n:]
            if df.shape[0] > 0:
                print(colors.blue(colors.bold(orm.__name__)))
                display(df)
