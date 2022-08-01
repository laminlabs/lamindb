import pandas as pd
from lndb_setup import settings

from ..dev import filepath_from_dobject
from ..dev.file import load_to_memory
from ..schema import core


class load:
    """Load data."""

    @classmethod
    def entity(cls, entity_name) -> pd.DataFrame:
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

    @classmethod
    def dobject(cls, dobject: core.dobject):
        """Load dobject into memory."""
        filepath = filepath_from_dobject(dobject)
        return load_to_memory(filepath)
