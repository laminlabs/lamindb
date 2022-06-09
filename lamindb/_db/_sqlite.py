import random
import string
from pathlib import Path

import pandas as pd
import sqlalchemy as sql

from .._settings import settings


def id(n_char: int = 6):
    base62 = string.digits + string.ascii_letters.swapcase()
    id = "".join(random.choice(base62) for i in range(n_char))
    return id


def id_file():
    return id(n_char=20)


def id_user():
    return id(n_char=3)


def get_database_file() -> Path:
    database_filename = str(settings.storage_root.stem).lower()  # type: ignore
    return settings.storage_root / f"{database_filename}.lndb"  # type: ignore


def get_engine():
    database_file = get_database_file()
    return sql.create_engine(f"sqlite:///{database_file}", future=True)


class insert_if_not_exists:
    """Insert data if it does not yet exist."""

    @classmethod
    def user(cls, user_name):
        df_user = db.load("user")

        if user_name in df_user.name.values:
            user_id = df_user.index[df_user.name == user_name][0]
            print(f"user {user_name} ({user_id}) already exists")
        else:
            user_id = db.insert.user(user_name)  # type: ignore
            print(f"added user {user_name} ({user_id})")

        return user_id


class insert:
    """Insert data."""

    @classmethod
    def user(cls, user_name):
        """User."""
        engine = get_engine()
        metadata = sql.MetaData()

        user = sql.Table(
            "user",
            metadata,
            sql.Column("id", sql.String, primary_key=True, default=id_user),
            autoload_with=engine,
        )

        with engine.begin() as conn:
            stmt = sql.insert(user).values(name=user_name)
            result = conn.execute(stmt)

        return result.inserted_primary_key[0]

    @classmethod
    def file(cls, name: str, *, source: str = None, source_name: str = None):
        """Data file with its origin."""
        source_id = source
        engine = get_engine()
        metadata = sql.MetaData()

        source_table = sql.Table(
            "file_source",
            metadata,
            autoload_with=engine,
        )

        file_table = sql.Table(  # primary key gen does not work with reflecting
            "file",
            metadata,
            sql.Column("id", sql.String, primary_key=True, default=id_file),
            autoload_with=engine,
        )

        if source_id is None:
            from nbproject import meta

            source_id = meta.id
            source_name = meta.title
            source_dependency = meta.dependency
            source_type = "nbproject"

            if source_name is None:
                raise RuntimeError(
                    "Can only ingest from notebook with title. Please set a title!"
                )
        else:
            source_dependency = None
            source_type = "other"

        from lamindb._configuration import user_id, user_name

        df_source = db.load("file_source")
        if source_id not in df_source.index:
            with engine.begin() as conn:
                stmt = sql.insert(source_table).values(
                    id=source_id,
                    name=source_name,
                    dependency=source_dependency,
                    type=source_type,
                    user=user_id,
                )
                conn.execute(stmt)
                print(
                    f"added source {source_name!r} ({source_id}) by user"
                    f" {user_name} ({user_id})"
                )

        with engine.begin() as conn:
            stmt = sql.insert(file_table).values(
                name=name,
                source=source_id,
            )
            result = conn.execute(stmt)
            file_id = result.inserted_primary_key[0]

        return file_id


class meta:
    """Manipulate the DB schema."""

    @classmethod
    def create(cls) -> None:
        """Create database with initial schema."""
        if get_database_file().exists():
            print("database already exists")
            return None

        # use the schema we just migrated to SQL and add a primary key
        metadata = sql.MetaData()

        # a user operating the database, e.g., ingesting data
        sql.Table(
            "user",
            metadata,
            sql.Column("id", sql.String, primary_key=True, default=id_user),
            sql.Column("name", sql.String),
        )

        # the entity that ingests the data file, the source of the data file
        # where the data file comes from
        # can be a notebook or a script/pipeline
        sql.Table(
            "file_source",
            metadata,
            sql.Column("id", sql.String, primary_key=True),  # this is an nbproject id
            sql.Column("name", sql.String),
            sql.Column("dependency", sql.String),
            sql.Column("type", sql.String),
            sql.Column("user", sql.String, sql.ForeignKey("user.id")),
        )

        # the data file
        sql.Table(
            "file",
            metadata,
            sql.Column("id", sql.String, primary_key=True, default=id_file),
            sql.Column("name", sql.String),
            sql.Column("source", sql.ForeignKey("file_source.id")),
        )

        engine = get_engine()
        metadata.create_all(bind=engine)

        print(f"created database at {get_database_file()}")


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
    12     >3e+21 (nbproject id)
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
    def insert_if_not_exists(cls):
        """Insert data if it does not exist."""
        return insert_if_not_exists

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
    def load(cls, entity_name) -> pd.DataFrame:
        """Load observations of entity as dataframe."""
        engine = get_engine()
        with engine.connect() as conn:
            df = pd.read_sql_table(entity_name, conn, index_col="id")
        return df
