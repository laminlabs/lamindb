from lamin_logger import logger
from lndb_cli._settings import (  # noqa
    load_or_create_instance_settings,
    load_or_create_user_settings,
)
from sqlmodel import SQLModel

from lamindb.dev.db import insert, insert_if_not_exists
from lamindb.dev.db._engine import get_engine
from lamindb.do import load


def setup_instance_db():
    """Setup database.

    Contains:
    - Database creation.
    - Sign-up and/or log-in.
    """
    instance_settings = load_or_create_instance_settings()
    user_settings = load_or_create_user_settings()
    if instance_settings.storage_dir is None:
        logger.warning("Instance is not configured. Call `lndb init` or `lndb load`.")
        return None
    instance_name = instance_settings.instance_name
    sqlite_file = instance_settings._sqlite_file
    from lamindb_schema import __version__

    if sqlite_file.exists():
        logger.info(f"Using instance: {sqlite_file}")
        schema_version = load("schema_version")
        if schema_version.index[-1] != __version__:
            result = input(
                "[Run this dialogue on the commmand line *outside* Jupyter]\nDid you"
                f" already migrate your db to schema version {__version__}? (y/n)"
            )
            if result == "y":
                insert.schema_version(__version__, user_settings.user_id)
            else:
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
