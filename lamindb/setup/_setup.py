from pathlib import Path
from typing import Union

from appdirs import AppDirs
from cloudpathlib import CloudPath
from sqlmodel import SQLModel

from lamindb.admin.db._engine import get_engine

from .._logger import logger
from ..dev._docs import doc_args
from ._settings import Settings, SettingsStore, _write, description, settings_file

DIRS = AppDirs("lamindb", "laminlabs")


def setup_storage_root(storage: Union[str, Path, CloudPath]) -> Union[Path, CloudPath]:
    if str(storage).startswith(("s3://", "gs://")):
        storage_root = CloudPath(storage)
    else:
        storage_root = Path(storage)
        if not storage_root.exists():
            storage_root.mkdir(parents=True)

    return storage_root


def setup_cache_root(
    settings: Settings,
) -> Union[Path, None]:
    if settings.cloud_storage:
        cache_root = Path(DIRS.user_cache_dir)
        if not cache_root.exists():
            cache_root.mkdir(parents=True)
    else:
        cache_root = None

    return cache_root


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


def setup_db(user_name):
    """Setup database.

    Contains:
    - Database creation.
    - Sign-up and/or log-in.
    """
    from lamindb.admin.db import insert_if_not_exists

    def create_db() -> None:
        """Create database with initial schema."""
        if settings()._db_file.exists():
            return None

        SQLModel.metadata.create_all(get_engine())

        logger.info(f"Created database {settings().db}.")

    create_db()

    user_id = insert_if_not_exists.user(user_name)

    return user_id


@doc_args(description.storage_root, description.user_name)
def setup_from_cli(
    *,
    storage: str,
    user: str,
) -> None:
    """Setup LaminDB. Alternative to using the CLI via `lamindb setup`.

    Args:
        storage: {}
        user: {}
    """
    settings = Settings()

    settings.user_name = user
    settings.storage_root = setup_storage_root(storage)
    settings.cache_root = setup_cache_root(settings)

    _write(settings)

    settings.user_id = setup_db(settings.user_name)

    _write(settings)  # type: ignore


def setup_from_store(store: SettingsStore) -> Settings:
    settings = Settings()

    settings.user_name = store.user_name
    settings.storage_root = setup_storage_root(store.storage_root)
    settings.cache_root = Path(store.cache_root) if store.cache_root != "null" else None
    settings.user_id = store.user_id

    return settings


def settings() -> Settings:
    """Return current settings."""
    if not settings_file.exists():
        logger.warning("Please setup lamindb via the CLI: lamindb setup")
        global Settings
        return Settings()
    else:
        settings_store = SettingsStore(_env_file=settings_file)
        settings = setup_from_store(settings_store)
        return settings
