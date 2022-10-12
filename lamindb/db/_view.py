from IPython.display import display

from .. import schema
from ..dev.db._select import select


def view(n: int = 10):
    """View the latest edits to the database.

    Args:
        n: display the latest n rows of the tables
    """
    modules = [
        i
        for i in schema.__dir__()
        if getattr(schema, i).__class__.__name__ == "module" and not i.startswith("_")
    ]

    # core tables always display first
    for module in ["core"] + [i for i in modules if i != "core"]:
        divider = "*" * 5
        print(f"{divider} {module} {divider}")
        schema_module = getattr(schema, module)
        tables = [
            i
            for i in schema_module.__dir__()
            if getattr(schema_module, i).__class__.__name__ == "SQLModelMetaclass"
        ]
        for entity in tables:
            df = getattr(select, entity)().df()
            if df.shape[0] > 0:
                print(entity)
                display(df.iloc[-n:])
