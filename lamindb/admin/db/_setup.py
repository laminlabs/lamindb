from sqlmodel import SQLModel

from . import get_database_file, get_engine


def setup() -> None:
    """Create database with initial schema."""
    if get_database_file().exists():
        print("database already exists")
        return None

    SQLModel.metadata.create_all(get_engine())

    print(f"created database at {get_database_file()}")
