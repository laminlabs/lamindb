from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath
from sqlmodel import SQLModel

from lamindb.admin.db import insert, insert_if_not_exists
from lamindb.admin.db._engine import get_engine
from lamindb.do import load

from .._logger import logger
from ..dev._docs import doc_args
from ._hub import sign_in_hub, sign_up_hub
from ._settings import (
    description,
    load_instance_settings,
    load_user_settings,
    setup_storage_dir,
    write_instance_settings,
    write_user_settings,
)


def setup_instance_db():
    """Setup database.

    Contains:
    - Database creation.
    - Sign-up and/or log-in.
    """
    instance_settings = load_instance_settings()
    user_settings = load_user_settings()
    instance_name = instance_settings.instance_name
    sqlite_file = instance_settings._sqlite_file
    from lamindb_schema import __version__

    if sqlite_file.exists():
        logger.info(f"Using instance: {sqlite_file}")
        schema_version = load("schema_version")
        if schema_version.index[-1] != __version__:
            raise RuntimeError(
                "\nEither migrate your instance db schema to version"
                f" {__version__}.\nOr install pip package lamindb_schema version"
                f" {schema_version.index[-1]}."
            )
    else:
        SQLModel.metadata.create_all(get_engine())
        instance_settings._update_cloud_sqlite_file()
        insert.schema_version(__version__, user_settings.user_id)
        logger.info(f"Created instance {instance_name}: {sqlite_file}")

    insert_if_not_exists.user(user_settings.user_email, user_settings.user_id)


def sign_up_first_time(email):
    user_settings = load_user_settings()
    user_settings.user_email = email
    write_user_settings(user_settings)
    secret = sign_up_hub(email)
    if secret is None:  # user already exists
        raise RuntimeError(
            "\nUser already exists! Please login instead: `lndb login`\n"
        )
    user_settings.user_secret = secret
    write_user_settings(user_settings)
    return None  # user needs to confirm email now


def log_in_user(
    *,
    email: Union[str, None] = None,
    secret: Union[str, None] = None,
):
    user_settings = load_user_settings()

    # user_email
    if email is None:
        if user_settings.user_email is None:
            raise RuntimeError(
                "No stored user email, please call: lndb login --email <your-email>"
            )
        else:
            email = user_settings.user_email
    else:
        user_settings.user_email = email
    write_user_settings(user_settings)

    # user_secret
    if secret is None:
        if user_settings.user_secret is None:
            raise RuntimeError(
                "No stored user secret, please call: lndb login --email <your-email>"
                " --email <your-secret>"
            )
        else:
            secret = user_settings.user_secret
    else:
        user_settings.user_secret = secret
    write_user_settings(user_settings)

    # user_id
    user_id = sign_in_hub(user_settings.user_email, user_settings.user_secret)
    user_settings.user_id = user_id
    write_user_settings(user_settings)


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
    instance_settings = load_instance_settings()
    user_settings = load_user_settings()
    if user_settings.user_id is None:
        if (
            user_settings.user_email is not None
            and user_settings.user_secret is not None  # noqa
        ):
            # complete user setup, this *only* happens after *sign_up_first_time*
            logger.info("Completing user sign up. Only happens once!")
            log_in_user(
                email=user_settings.user_email, secret=user_settings.user_secret
            )
            user_settings = load_user_settings()  # need to reload, here, to get user_id
        else:
            raise RuntimeError("Login user: lndb login --email")
    write_user_settings(user_settings)

    # setup storage
    if storage is None:
        if instance_settings.storage_dir is None:
            raise RuntimeError(
                "No storage in .env, please call: lndb init --storage <location>"
            )
        else:
            storage = instance_settings.storage_dir
    else:
        instance_settings.storage_dir = setup_storage_dir(storage)
    write_instance_settings(instance_settings)

    # setup _config
    instance_settings._dbconfig = dbconfig
    if dbconfig != "sqlite":
        raise NotImplementedError()
    write_instance_settings(instance_settings)

    # setup instance db
    setup_instance_db()
    return None
