from __future__ import annotations

import builtins
import importlib
import inspect

from lamin_utils import colors, logger
from lamindb_setup import settings
from lamindb_setup._init_instance import get_schema_module_name
from lnschema_core import Record

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)


def view(
    n: int = 7, schema: str | None = None, registries: list[str] | None = None
) -> None:
    """View latest metadata state.

    Args:
        n: Display the last `n` rows of a registry.
        schema: Schema module to view. Default's to
            `None` and displays all schema modules.
        registries: List of Record names. Defaults to
            `None` and lists all registries.

    Examples:
        >>> ln.view()
    """
    if is_run_from_ipython:
        from IPython.display import display as show
    else:
        show = logger.print

    if schema is not None:
        schema_names = [schema]
    else:
        schema_names = ["core"] + list(settings.instance.schema)

    for schema_name in schema_names:
        schema_module = importlib.import_module(get_schema_module_name(schema_name))

        all_registries = {
            orm
            for orm in schema_module.__dict__.values()
            if inspect.isclass(orm)
            and issubclass(orm, Record)
            and orm.__name__ != "Record"
        }
        if registries is not None:
            filtered_registries = {
                orm for orm in all_registries if orm.__name__ in registries
            }
        else:
            filtered_registries = all_registries
        if len(schema_names) > 1:
            section = f"* module: {colors.green(colors.bold(schema_name))} *"
            section_no_color = f"* module: {schema_name} *"
            logger.print("*" * len(section_no_color))
            logger.print(section)
            logger.print("*" * len(section_no_color))
        for orm in sorted(filtered_registries, key=lambda x: x.__name__):
            if hasattr(orm, "updated_at"):
                df = orm.filter().order_by("-updated_at")[:n].df()
            else:
                # need to adjust in the future
                df = orm.df().iloc[-n:]
            if df.shape[0] > 0:
                logger.print(colors.blue(colors.bold(orm.__name__)))
                show(df)
