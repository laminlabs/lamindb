import random
import sqlite3
import string
from pathlib import Path
from subprocess import call

import pandas as pd
from sqlalchemy import (
    Column,
    ForeignKey,
    MetaData,
    String,
    Table,
    create_engine,
    insert,
    update,
)


def uid(n_char=6):
    base62 = string.digits + string.ascii_letters.swapcase()
    id = "".join(random.choice(base62) for i in range(n_char))
    return id


def ingest(entity_name: str, df: pd.DataFrame, *, len_id: int = 6) -> None:
    """Ingest data for an entity into the database.

    Args:
        entity_name: Name of the entity.
        df: DataFrame storing data about observations of the entity.
        len_id: Length of the base62 ID uniquely characterizing the entity.

    The following table shows the number of entries available given `len_id`.

    ====== =========
    len_id n_entries
    ====== =========
    1      >6e+01
    2      >4e+03
    3      >2e+05
    4      >1e+07
    5      >9e+08
    6      >6e+10
    7      >4e+12
    8      >2e+14
    9      >1e+16
    ====== =========
    """
    memory = create_engine("sqlite:///:memory:", future=True)

    if df.duplicated().any():
        raise ValueError("df contains duplicate rows, pass df.drop_duplicates()")

    _df = df.reset_index()  # we need the index for relations

    # this also writes the data but we ignore this for now
    with memory.begin() as conn:
        _df.to_sql(entity_name, conn, index=False)

    def _uid():
        return uid(n_char=len_id)

    # use the schema we just migrated to SQL and add a primary key
    metadata = MetaData()
    table = Table(
        entity_name,
        metadata,
        Column("id", String, primary_key=True, default=_uid),
        autoload_with=memory,
    )

    engine = create_engine("sqlite:///graham.db", future=True)
    metadata.create_all(bind=engine)

    with engine.connect() as conn:
        for row in _df.iterrows():
            conn.execute(insert(table), row[1].to_dict())
        conn.commit()


def create(entity_name: str, df: pd.DataFrame, *, len_id: int = 6, seed=0) -> None:
    """Create the database starting with data for one entity.

    Args:
        entity_name: Name of the entity.
        df: DataFrame storing data about observations of the entity.
        len_id: Length of the base62 ID uniquely characterizing the entity.
        seed: Seed for generating pseudo-random hash-based IDs.

    The following table shows the number of entries available given `len_id`.

    ====== =========
    len_id n_entries
    ====== =========
    1      >6e+01
    2      >4e+03
    3      >2e+05
    4      >1e+07
    5      >9e+08
    6      >6e+10
    7      >4e+12
    8      >2e+14
    9      >1e+16
    ====== =========
    """
    if Path("graham.db").exists():
        call("rm graham.db", shell=True)
    random.seed(seed)
    return ingest(entity_name, df, len_id=len_id)


def load(entity_name, drop_index=True) -> pd.DataFrame:
    """Load observations of entity as dataframe."""
    engine = create_engine("sqlite:///graham.db", future=True)
    with engine.begin() as conn:
        df = pd.read_sql_table(entity_name, conn)
    if drop_index:
        df = df.drop(columns=["index"])
    return df.set_index("id")


def entities():
    """Return all entities in the db."""
    metadata = MetaData()
    engine = create_engine("sqlite:///graham.db", future=True)
    metadata.reflect(bind=engine)
    table_names = [table.name for table in metadata.sorted_tables]
    return table_names


def link(reference, target, df) -> None:
    """Link two entities.

    Args:
        reference: Entity name of reference entity.
        target: Entity name of target entity.
        df: DataFrame with linkage information.
    """
    df_linkage = df
    df_reference = load(reference)
    if df_reference.shape[0] != df_linkage.shape[0]:
        raise ValueError("df for linkage needs to have same length as reference table")
    df_target = load(target)
    df_target[target] = df_target.index
    df_linkage = df_linkage.merge(df_target, validate="many_to_one")

    metadata = MetaData()
    engine = create_engine("sqlite:///graham.db", future=True)
    metadata.create_all(bind=engine)

    table = Table(
        reference,
        metadata,
        Column(target, String, ForeignKey("species.id")),
        autoload_with=engine,
    )

    metadata.create_all(bind=engine)

    with sqlite3.connect("graham.db") as conn:
        cur = conn.cursor()
        stmt = f"ALTER TABLE {reference} ADD COLUMN {target} REFERENCES {target}(id)"
        cur.execute(stmt)

    with engine.connect() as conn:
        for row in df_linkage.iterrows():
            conn.execute(
                update(table).where(table.c.index == row[0]), {target: row[1][target]}
            )
        conn.commit()
