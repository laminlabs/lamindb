from sqlmodel import SQLModel


def create_db() -> None:
    """Create database with initial schema."""
    from lamindb.admin.db._engine import get_engine

    from . import settings

    if settings._db_file.exists():
        print("database already exists")
        return None

    SQLModel.metadata.create_all(get_engine())

    print(f"created database at {settings.db}")


def setup_db(settings):
    """Register user in database."""
    from lamindb.admin.db import insert_if_not_exists

    create_db()

    user_id = insert_if_not_exists.user(settings["user_name"])

    settings["user_id"] = user_id
