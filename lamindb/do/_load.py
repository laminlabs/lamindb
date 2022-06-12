import pandas as pd

from ._admin.db import get_engine


def load(entity_name) -> pd.DataFrame:
    """Load observations of entity as dataframe."""
    engine = get_engine()
    with engine.connect() as conn:
        df = pd.read_sql_table(entity_name, conn, index_col="id")
    return df
