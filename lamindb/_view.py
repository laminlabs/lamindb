from __future__ import annotations

import builtins
import importlib
import inspect
from typing import TYPE_CHECKING

from lamin_utils import colors, logger
from lamindb_setup import settings
from lamindb_setup._init_instance import get_schema_module_name

from lamindb.models import Feature, FeatureValue, ParamValue, SQLRecord

from .models.feature import serialize_pandas_dtype

if TYPE_CHECKING:
    import pandas as pd

is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)


def display_df_with_descriptions(
    df: pd.DataFrame, descriptions: dict[str, str] | None = None
):
    from IPython.display import HTML, display

    if descriptions is None:
        display(df)
        return None

    # Start building HTML table
    html = '<table class="dataframe">'

    # Create header with title and description rows
    html += "<thead>"

    # Column names row
    html += "<tr>"
    html += '<th class="header-title index-header"></th>'  # Index header
    for col in df.columns:
        html += f'<th class="header-title">{col}</th>'
    html += "</tr>"

    # Descriptions row
    html += "<tr>"
    html += f'<th class="header-desc index-header">{df.index.name or ""}</th>'  # Index column
    for col in df.columns:
        desc = descriptions.get(col, "")
        html += f'<th class="header-desc">{desc}</th>'
    html += "</tr>"

    html += "</thead>"

    # Add body rows
    html += "<tbody>"
    for idx, row in df.iterrows():
        html += "<tr>"
        html += f'<th class="row-index">{idx}</th>'  # Index value
        for col in df.columns:
            html += f"<td>{row[col]}</td>"
        html += "</tr>"
    html += "</tbody>"
    html += "</table>"

    # Add CSS styles
    styled_html = f"""
    <style>
        .dataframe {{
            border-collapse: collapse;
            margin: 10px 0;
        }}
        .dataframe th, .dataframe td {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }}
        .header-title {{
            font-weight: bold;
        }}
        .header-desc {{
            color: #666;
            font-weight: normal;
        }}
        .row-index {{
            font-weight: bold;
        }}
        .index-header {{
            font-weight: bold;
        }}
    </style>
    {html}
    """
    return display(HTML(styled_html))


def view(
    df: pd.DataFrame | None = None,
    limit: int = 7,
    modules: str | None = None,
    registries: list[str] | None = None,
) -> None:
    """View metadata.

    Args:
        df: A DataFrame to display.
        limit: Display the latest `n` records
        modules: schema module to view. Default's to
            `None` and displays all registry modules.
        registries: List of SQLRecord names. Defaults to
            `None` and lists all registries.

    Examples:
        >>> ln.view()
    """
    if df is not None:
        descriptions = {
            col_name: serialize_pandas_dtype(dtype)
            for col_name, dtype in df.dtypes.to_dict().items()
        }
        feature_dtypes = dict(Feature.objects.values_list("name", "dtype"))
        descriptions.update(feature_dtypes)
        display_df_with_descriptions(df, descriptions)
        return None

    if is_run_from_ipython:
        from IPython.display import display as show
    else:
        show = logger.print

    if modules is not None:
        module_names = [modules]
    else:
        module_names = ["core"] + list(settings.instance.modules)

    for module_name in module_names:
        schema_module = importlib.import_module(get_schema_module_name(module_name))
        # the below is necessary because a schema module might not have been
        # explicitly accessed
        importlib.reload(schema_module)

        all_registries = {
            registry
            for registry in schema_module.__dict__.values()
            if inspect.isclass(registry)
            and issubclass(registry, SQLRecord)
            and registry is not SQLRecord
        }
        if module_name == "core":
            all_registries.update({FeatureValue, ParamValue})
        if registries is not None:
            filtered_registries = {
                registry
                for registry in all_registries
                if registry.__name__ in registries
            }
        else:
            filtered_registries = all_registries
        if len(module_names) > 1:
            section = f"* module: {colors.green(colors.bold(module_name))} *"
            section_no_color = f"* module: {module_name} *"
            logger.print("*" * len(section_no_color))
            logger.print(section)
            logger.print("*" * len(section_no_color))
        for registry in sorted(filtered_registries, key=lambda x: x.__name__):
            df = registry.df(limit=limit)
            if df.shape[0] > 0:
                logger.print(colors.blue(colors.bold(registry.__name__)))
                show(df)
