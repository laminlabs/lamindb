from pathlib import Path
from typing import Union

from appdirs import AppDirs
from cloudpathlib import CloudPath
from sqlmodel import SQLModel

from lamindb.admin.db import insert_if_not_exists
from lamindb.admin.db._engine import get_engine

from .._logger import logger
from ..dev._docs import doc_args
from ._hub import sign_in_hub, sign_up_hub
from ._settings import (
    Settings,
    description,
    instance_from_storage,
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


def setup_instance_db():
    """Setup database.

    Contains:
    - Database creation.
    - Sign-up and/or log-in.
    """
    settings = load_settings()
    sqlite_file = settings._sqlite_file
    if sqlite_file.exists():
        logger.info(f"Using instance DB: {sqlite_file}")
    else:
        SQLModel.metadata.create_all(get_engine())
        settings._update_cloud_sqlite_file()
        logger.info(f"Created instance DB: {sqlite_file}")

    insert_if_not_exists.user(settings.user_email, settings.user_id)


def setup_user_and_secret(user, secret):
    # retrieves user email & secret from stored settings
    # also signs up new user to hub
    settings = load_settings()
    # user
    if user is None:
        if settings.user_email is None:
            raise RuntimeError(
                "No stored user email, please call: lndb setup --user <your-email>"
            )
        else:
            user = settings.user_email

    # secret
    if secret is None:
        if settings.user_secret is None:  # typically, secret is stored
            secret = sign_up_hub(user)
            if secret is None:  # user already exists
                raise RuntimeError(
                    "Please pass secret: lndb setup --secret <your-secret>"
                    "If you don't have it anymore, recover it via email: "
                    "https://hub.lamin.ai"
                )
            else:  # successful sign up of new user
                settings.user_secret = secret
                write_settings(settings)  # store user_email & secret
        else:
            secret = settings.user_secret
    else:
        settings.user_secret = secret
        write_settings(settings)  # update stored secret


@doc_args(
    description.user_email,
    description.user_secret,
    description.storage_dir,
    description.instance_name,
    description.db,
)
def setup_from_cli(
    *,
    user: str,
    secret: Union[str, None] = None,
    storage: Union[str, Path, CloudPath],
    instance: Union[str, None] = None,
    db: str = "sqlite",
) -> None:
    """Setup LaminDB.

    Args:
        user: {}
        secret: {}
        instance: {}
        storage: {}
        db: {}
    """
    # the following ensures that user & secret are stored in the settings .env
    setup_user_and_secret(user, secret)

    # settings.user_email & settings.user_secret are set
    settings = load_settings()

    # continue with signing into hub
    user_id = sign_in_hub(settings.user_email, settings.user_secret)
    settings.user_id = user_id
    write_settings(settings)

    # setup storage
    if storage is None:
        if settings.storage_dir is None:
            raise RuntimeError(
                "No storage in .env, please call: lndb setup --storage <location>"
            )
        else:
            storage = settings.storage_dir
    else:
        settings.storage_dir = setup_storage_dir(storage)
    write_settings(settings)

    # setup instance_name
    if instance is None:
        if settings.instance_name is None:
            if db == "sqlite":
                instance = instance_from_storage(storage)
                logger.info(f"Inferred instance from storage: {instance}")
            else:
                raise RuntimeError("Provide an instance name!")
        else:
            instance = settings.instance_name
    else:
        if db == "sqlite":
            is_consistent = instance == instance_from_storage(storage)
            if not is_consistent:
                raise RuntimeError(
                    "Your instance name is determined by your storage location!\n"
                    f"It is auto-determined as: {instance_from_storage(storage)}\n"
                    "Do *not* pass `--instance` in `lndb set`."
                )

    settings.instance_name = instance
    write_settings(settings)

    # setup cache_dir
    settings.cache_dir = setup_cache_dir(settings)
    write_settings(settings)

    # setup instance db
    setup_instance_db()
    return None
