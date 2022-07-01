from sqlmodel import SQLModel

from lamindb.admin.db._engine import get_engine

from .._logger import logger
from . import settings


def create_db() -> None:
    """Create database with initial schema."""
    if settings()._db_file.exists():
        return None

    SQLModel.metadata.create_all(get_engine())

    logger.info(f"Created database at {settings().db}.")


def setup_db(user_name):
    """Register user in database."""
    from lamindb.admin.db import insert_if_not_exists

    create_db()

    user_id = insert_if_not_exists.user(user_name)

    return user_id
