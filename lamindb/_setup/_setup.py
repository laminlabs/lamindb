from pathlib import Path
from typing import Union

from appdirs import AppDirs
from cloudpathlib import CloudPath
from sqlmodel import SQLModel

from lamindb.admin.db._engine import get_engine

from .._logger import logger
from ..dev._docs import doc_args
from ._settings import (
    Settings,
    description,
    load_settings,
    setup_storage_dir,
    write_settings,
)

DIRS = AppDirs("lamindb", "laminlabs")


def setup_cache_dir(
    settings: Settings,
) -> Union[Path, None]:
    if settings.cloud_storage:
        cache_dir = Path(DIRS.user_cache_dir)
        if not cache_dir.exists():
            cache_dir.mkdir(parents=True)
    else:
        cache_dir = None

    return cache_dir


def setup_db(user_email, user_id):
    """Setup database.

    Contains:
    - Database creation.
    - Sign-up and/or log-in.
    """
    from lamindb.admin.db import insert_if_not_exists

    def create_db() -> None:
        """Create database with initial schema."""
        settings = load_settings()
        sqlite_file = settings._sqlite_file
        if sqlite_file.exists():
            logger.info(f"Using DB instance {settings.instance_name} at {sqlite_file}")
            return None

        SQLModel.metadata.create_all(get_engine())

        logger.info(f"Created DB instance {settings.instance_name} at {sqlite_file}.")
        settings._update_cloud_sqlite_file()

    create_db()

    user_id = insert_if_not_exists.user(user_email, user_id)

    return user_id


@doc_args(
    description.user_email,
    description.user_secret,
    description.storage_dir,
    description.instance_name,
)
def setup_from_cli(
    *,
    user: str,
    secret: Union[str, None] = None,
    storage: Union[str, Path, CloudPath],
    instance: Union[str, None] = None,
) -> None:
    """Setup LaminDB. Alternative to using the CLI via `lamindb setup`.

    Args:
        user: {}
        secret: {}
        storage: {}
        instance: {}
    """
    settings = load_settings()
    if user is None:
        if settings.user_email is None:
            user = input(f"Please paste your {description.user_email}: ")
            settings.user_email = user
        else:
            user = settings.user_email

    if secret is None:
        if settings.user_secret is None:
            user = input(f"Please paste your {description.user_secret}: ")
            settings.user_secret = secret
        else:
            secret = settings.user_secret

    settings.storage_dir = setup_storage_dir(storage)

    # setup instance
    if instance is None:
        storage = str(storage)
        if storage.startswith(("s3://", "gs://")):
            instance = storage.replace("s3://", "")
        else:
            instance = str(Path(storage).stem)
    settings.instance_name = instance

    # setup cache_dir
    settings.cache_dir = setup_cache_dir(settings)

    write_settings(settings)

    if secret is None:  # sign up
        settings.user_secret = setup_db(settings.user_email, secret)
        write_settings(settings)  # type: ignore
    else:  # log in
        settings.user_id = setup_db(settings.user_email, secret)
        write_settings(settings)  # type: ignore

    return None
