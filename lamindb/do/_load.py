import pandas as pd
from lndb_cli import load_or_create_instance_settings


def load(entity_name) -> pd.DataFrame:
    """Load observations of entity as dataframe."""
    settings = load_or_create_instance_settings()
    engine = settings.db_engine()
    with engine.connect() as conn:
        df = pd.read_sql_table(entity_name, conn)
        if "id" in df.columns:
            if "v" in df.columns:
                df = df.set_index(["id", "v"])
            else:
                df = df.set_index("id")
    return df
