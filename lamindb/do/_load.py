import pandas as pd
from lndb_setup import settings


def load(entity_name) -> pd.DataFrame:
    """Load observations of entity as dataframe."""
    engine = settings.instance.db_engine()
    with engine.connect() as conn:
        df = pd.read_sql_table(entity_name, conn)
        if "id" in df.columns:
            if "v" in df.columns:
                df = df.set_index(["id", "v"])
            else:
                df = df.set_index("id")
    return df
