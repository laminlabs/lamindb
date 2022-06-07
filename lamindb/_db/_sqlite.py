import random
import string
from pathlib import Path

import pandas as pd
import sqlalchemy as sql
from sqlalchemy import Column, String

from .._settings import settings


def uid(n_char=6):
    base62 = string.digits + string.ascii_letters.swapcase()
    id = "".join(random.choice(base62) for i in range(n_char))
    return id


def uid_file():
    return uid(n_char=20)


def uid_user():
    return uid(n_char=3)


def get_database_file() -> Path:
    database_filename = str(settings.storage_root.stem).lower()  # type: ignore
    return settings.storage_root / f"{database_filename}.lndb"  # type: ignore


def get_engine():
    database_file = get_database_file()
    return sql.create_engine(f"sqlite:///{database_file}", future=True)


class insert:
    """Insert data."""

    @classmethod
    def file(cls, name, source):
        """Data file."""
        engine = get_engine()
        metadata = sql.MetaData()

        file = sql.Table(  # primary key gen does not work with reflecting
            "file",
            metadata,
            Column("id", String, primary_key=True, default=uid_file),
            autoload_with=engine,
        )

        notebook = sql.Table(
            "source",
            metadata,
            autoload_with=engine,
        )

        with engine.begin() as conn:
            stmt = sql.insert(file).values(
                name=name,
                source=source,
            )
            result = conn.execute(stmt)

            from lamindb._configuration import user_id

            sql.insert(notebook).values(
                id=source,
                user=user_id,
            )
            conn.execute(stmt)

        return result.inserted_primary_key[0]

    @classmethod
    def user(cls):
        """User."""
        engine = get_engine()
        metadata = sql.MetaData()

        user = sql.Table(
            "user",
            metadata,
            Column("id", String, primary_key=True, default=uid_user),
            autoload_with=engine,
        )

        from lamindb._configuration import user_name

        with engine.begin() as conn:
            stmt = sql.insert(user).values(name=user_name)
            result = conn.execute(stmt)

        return result.inserted_primary_key[0], user_name


class meta:
    """Manipulate the DB schema."""

    @classmethod
    def create(cls) -> None:
        """Create database with initial schema."""
        if get_database_file().exists():
            print("database already exists")
            # add a check for whether the user already exists!
            user_id, user_name = db.insert.user()  # type: ignore
            print(f"adding user {user_id} ({user_name})")
            return None

        # use the schema we just migrated to SQL and add a primary key
        metadata = sql.MetaData()

        # the data file
        sql.Table(
            "file",
            metadata,
            Column("id", String, primary_key=True, default=uid_file),
            Column("name", String),
            Column("source", sql.ForeignKey("source.id")),
        )

        # the entity that ingests the data file, the source of the data file
        # where the data file comes from
        # can be a notebook or a script/pipeline
        sql.Table(
            "source",
            metadata,
            Column("id", String, primary_key=True),  # this is an nbproject uid
            Column("name", String),
            Column("user", String, sql.ForeignKey("user.id")),
        )

        # a user operating on the database, e.g., ingesting data
        sql.Table(
            "user",
            metadata,
            Column("id", String, primary_key=True, default=uid_user),
            Column("name", String),  # can be anything, e.g. a URL
        )

        engine = get_engine()
        metadata.create_all(bind=engine)

        user_id, user_name = db.insert.user()  # type: ignore

        print(
            f"created database at {get_database_file()} by user {user_id} ({user_name})"
        )
        return user_id


class db:
    """Interact with local database & storage.

    We work with base62 primary keys:

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
    12     >3e+21 (nbproject uid)
    20     >7e+35 (~UUID)
    ====== =========
    """

    @classmethod
    @property
    def insert(cls):
        """Insert data."""
        return insert

    @classmethod
    @property
    def meta(cls):
        """Change the schema, create database."""
        return meta

    @classmethod
    @property
    def entities(cls):
        """Return all entities in the db."""
        metadata = sql.MetaData()
        engine = get_engine()
        metadata.reflect(bind=engine)
        table_names = [table.name for table in metadata.sorted_tables]
        return table_names

    @classmethod
    def load(entity_name, drop_index=True) -> pd.DataFrame:
        """Load observations of entity as dataframe."""
        engine = get_engine()
        with engine.connect() as conn:
            df = pd.read_sql_table(entity_name, conn)
        if drop_index:
            df = df.drop(columns=["index"])
        return df.set_index("id")
