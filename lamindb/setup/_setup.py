import base64
from pathlib import Path
from typing import Union
from urllib.request import urlretrieve

from appdirs import AppDirs
from cloudpathlib import CloudPath
from sqlmodel import SQLModel

from lamindb.admin.db._engine import get_engine

from .._logger import logger
from ..dev import id
from ..dev._docs import doc_args
from ._settings import (
    Settings,
    description,
    load_settings,
    setup_storage_dir,
    write_settings,
)
from ._settings_store import Connector

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


def setup_notion(notion: str = None):
    NOTION_HELP = (
        "Notion integration token (currently dysfunctional, see"
        " https://lamin.ai/notes/2022/explore-notion)"
    )

    if notion is None:
        notion = input(
            f"Please paste your internal {NOTION_HELP}"
            " (https://notion.so/my-integrations): "
        )

    return dict(NOTION_API_KEY=notion)


def setup_db(user_email, secret=None):
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

    from supabase import create_client

    connector_file, _ = urlretrieve(
        "https://lamin-site-assets.s3.amazonaws.com/connector.env"
    )
    connector = Connector(_env_file=connector_file)

    supabase = create_client(connector.url, connector.key)

    if secret is None:
        secret = id.id_secret()
        supabase.auth.sign_up(email=user_email, password=secret)
        logger.info(
            f"Generated login secret: {secret}.\n"
            "Please confirm the sign-up email and then repeat the login call.\n"
            "Your secret has been stored and will persist until a new installation."
        )
        return secret
    else:
        session = supabase.auth.sign_in(email=user_email, password=secret)

    uuid_b64 = base64.urlsafe_b64encode(session.user.id.bytes).decode("ascii")
    user_id = uuid_b64[:8].replace("_", "0").replace("-", "1")

    # data = supabase.table("usermeta").insert(
    #     {"id": session.user.id.hex, "lnid": user_id}
    # ).execute()
    # assert len(data.data) > 0

    supabase.auth.sign_out()

    user_id = insert_if_not_exists.user(user_email, user_id)

    return user_id


@doc_args(description.storage_dir, description.user_email, description.instance_name)
def setup_from_cli(
    *,
    storage: Union[str, Path, CloudPath],
    user: str,
    secret: Union[str, None] = None,
    instance: Union[str, None] = None,
) -> None:
    """Setup LaminDB. Alternative to using the CLI via `lamindb setup`.

    Args:
        storage: {}
        user: {}
        secret: Generated secret for login.
        instance: {}
    """
    if secret is None:
        settings = load_settings()
        if settings.secret is not None:
            secret = settings.secret

    settings = Settings(secret=secret)

    # setup user & storage
    settings.user_email = user
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
        settings.secret = setup_db(settings.user_email, secret)
        write_settings(settings)  # type: ignore
    else:  # log in
        settings.user_id = setup_db(settings.user_email, secret)
        write_settings(settings)  # type: ignore

    return None
