import random
import string
from pathlib import Path

import sqlalchemy as sql
from sqlalchemy import Column, String

from .._settings import settings


def uid(n_char=6):
    base62 = string.digits + string.ascii_letters.swapcase()
    id = "".join(random.choice(base62) for i in range(n_char))
    return id


def uid_file():
    return uid(n_char=20)


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
            "notebook",
            metadata,
            autoload_with=engine,
        )

        with engine.begin() as conn:
            stmt = sql.insert(file).values(
                name=name,
                source=source,
            )
            result = conn.execute(stmt)
            sql.insert(notebook).values(
                id=source,
            )
            conn.execute(stmt)

        return result.inserted_primary_key[0]


class meta:
    """Manipulate the DB schema."""

    @classmethod
    def create(cls, seed: int = 0) -> None:
        """Create database with initial schema."""
        if get_database_file().exists():
            print("database already exists, create has no effect")
            return None

        # use the schema we just migrated to SQL and add a primary key
        metadata = sql.MetaData()

        sql.Table(
            "file",
            metadata,
            Column("id", String, primary_key=True, default=uid_file),
            Column("name", String),
            Column("source", String),  # can be anything, e.g. a notebook
        )

        sql.Table(
            "notebook",
            metadata,
            Column("id", String, primary_key=True),  # this is an nbproject uid
            Column("source", String),  # can be anything, e.g. a URL
        )

        engine = get_engine()
        metadata.create_all(bind=engine)

        print(f"created database at {get_database_file()}")
        return None


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
