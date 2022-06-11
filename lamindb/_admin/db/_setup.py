from typing import Optional

from sqlmodel import Field, SQLModel

from . import get_database_file, get_engine
from .id import id_file, id_user


def setup() -> None:
    """Create database with initial schema."""
    if get_database_file().exists():
        print("database already exists")
        return None

    # a user operating the database, e.g., ingesting data
    class user(SQLModel, table=True):  # type: ignore
        id: Optional[str] = Field(default=id_user, primary_key=True)
        name: str

    # the process that ingests the data file, the source of the data file
    # where the data file comes from
    # can be a notebook or a script/pipeline
    class file_source(SQLModel, table=True):  # type: ignore
        id: str = Field(default=None, primary_key=True)
        name: Optional[str]
        dependency: Optional[str]
        type: str
        user: str = Field(foreign_key="user.id")

    # the data file
    class file(SQLModel, table=True):  # type: ignore
        id: Optional[str] = Field(default=id_file, primary_key=True)
        name: str
        source: str = Field(foreign_key="file_source.id")

    SQLModel.metadata.create_all(get_engine())

    print(f"created database at {get_database_file()}")
