import sqlalchemy as sql

from . import get_database_file, get_engine
from .id import id_file, id_user


def setup() -> None:
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
