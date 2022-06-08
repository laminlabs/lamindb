import random
import string

import pandas as pd
import sqlalchemy as sql


def get_engine():
    return sql.create_engine("sqlite:///test.db", future=True)


def uid(n_char=4):
    base62 = string.digits + string.ascii_letters.swapcase()
    id = "".join(random.choice(base62) for i in range(n_char))
    return id


def create_with_sqlcore():

    metadata = sql.MetaData()

    source = sql.Table(
        "source",
        metadata,
        sql.Column("id", sql.String, primary_key=True, default=uid),
        sql.Column("name", sql.String),
    )

    file = sql.Table(
        "file",
        metadata,
        sql.Column("id", sql.String, primary_key=True, default=uid),
        sql.Column("name", sql.String),
        sql.Column("source", sql.ForeignKey("source.id")),
    )

    engine = get_engine()
    metadata.create_all(bind=engine)

    # insert data
    with engine.begin() as conn:
        stmt = sql.insert(source).values(
            name="Ingest Alpert19",
        )
        result = conn.execute(stmt)
        stmt = sql.insert(file).values(
            name="test_file.csv",
            source=result.inserted_primary_key[0],
        )
        conn.execute(stmt)


def load_table_pandas(table_name):
    engine = get_engine()
    with engine.connect() as conn:
        df = pd.read_sql_table(table_name, conn, index_col="id")
    return df


def test_sqlcore_to_pandas():

    create_with_sqlcore()
    print(load_table_pandas("file"))
    print(load_table_pandas("source"))


if __name__ == "__main__":
    test_sqlcore_to_pandas()
