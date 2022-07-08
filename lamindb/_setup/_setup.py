from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath
from sqlmodel import SQLModel

from lamindb.admin.db import insert_if_not_exists
from lamindb.admin.db._engine import get_engine

from .._logger import logger
from ..dev._docs import doc_args
from ._hub import sign_in_hub, sign_up_hub
from ._settings import description, load_settings, setup_storage_dir, write_settings


def setup_instance_db():
    """Setup database.

    Contains:
    - Database creation.
    - Sign-up and/or log-in.
    """
    settings = load_settings()
    instance_name = settings.instance_name
    sqlite_file = settings._sqlite_file
    if sqlite_file.exists():
        logger.info(f"Using instance: {sqlite_file}")
    else:
        SQLModel.metadata.create_all(get_engine())
        settings._update_cloud_sqlite_file()
        logger.info(f"Created instance {instance_name}: {sqlite_file}")

    insert_if_not_exists.user(settings.user_email, settings.user_id)


def sign_up_first_time(email):
    settings = load_settings()
    settings.user_email = email
    write_settings(settings)
    secret = sign_up_hub(email)
    if secret is None:  # user already exists
        raise RuntimeError(
            "\nUser already exists! Please login instead: `lndb login`\n"
        )
    settings.user_secret = secret
    write_settings(settings)
    return None  # user needs to confirm email now


def log_in_user(
    *,
    email: Union[str, None] = None,
    secret: Union[str, None] = None,
):
    settings = load_settings()

    # user_email
    if email is None:
        if settings.user_email is None:
            raise RuntimeError(
                "No stored user email, please call: lndb login --email <your-email>"
            )
        else:
            email = settings.user_email
    else:
        settings.user_email = email
    write_settings(settings)

    # user_secret
    if secret is None:
        if settings.user_secret is None:
            raise RuntimeError(
                "No stored user secret, please call: lndb login --email <your-email>"
                " --email <your-secret>"
            )
        else:
            secret = settings.user_secret
    else:
        settings.user_secret = secret
    write_settings(settings)

    # user_id
    user_id = sign_in_hub(settings.user_email, settings.user_secret)
    settings.user_id = user_id
    write_settings(settings)


@doc_args(
    description.storage_dir,
    description._dbconfig,
)
def setup_instance(
    *,
    storage: Union[str, Path, CloudPath],
    dbconfig: str = "sqlite",
) -> None:
    """Setup LaminDB.

    Args:
        storage: {}
        dbconfig: {}
    """
    # settings.user_email & settings.user_secret are set
    settings = load_settings()
    if settings.user_id is None:
        if settings.user_email is not None and settings.user_secret is not None:
            # complete user setup, this *only* happens after *sign_up_first_time*
            logger.info("Completing user sign up. Only happens once!")
            log_in_user(email=settings.user_email, secret=settings.user_secret)
            settings = load_settings()  # need to reload, here, to get user_id
        else:
            raise RuntimeError("Configure user: lndb user --email")

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

    # setup _config
    settings._dbconfig = dbconfig
    if dbconfig != "sqlite":
        write_settings(settings)
        raise NotImplementedError()

    # setup instance db
    setup_instance_db()
    return None
