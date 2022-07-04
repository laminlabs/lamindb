from dataclasses import dataclass
from pathlib import Path
from typing import Union, get_type_hints

from cloudpathlib import CloudPath, S3Client

from .._logger import logger
from ._settings_store import SettingsStore

root_dir = Path(__file__).parent.resolve()
settings_file = root_dir / ".env"


# A mere tool for quick access to the docstrings above
# I thought I had it work to read from the docstrings above, but doesn't seem so
class description:
    storage_dir = """Storage root. Either local dir, ``s3://bucket_name`` or ``gs://bucket_name``."""  # noqa
    cache_dir = """Cache root, a local directory to cache cloud files."""
    user_name = """User name. Consider using the GitHub username."""
    user_id = """User name. Consider using the GitHub username."""
    instance_name = """Name of LaminDB instance, which corresponds to exactly one backend database."""  # noqa


def storage_filepath(filekey: Union[Path, CloudPath, str]) -> Union[Path, CloudPath]:
    """Cloud or local filepath from filekey."""
    settings = load_settings()
    if settings.cloud_storage:
        client = S3Client(local_cache_dir=settings.cache_dir)
        return client.CloudPath(settings.storage_dir / filekey)
    else:
        return settings.storage_dir / filekey


def local(filepath: Union[Path, CloudPath]) -> Path:
    """Local (cache) filepath from filepath."""
    if load_settings().cloud_storage:
        filepath = filepath.fspath  # type: ignore  # mypy misses CloudPath
    return filepath


def local_filepath(filekey: Union[Path, CloudPath, str]) -> Path:
    """Local (cache) filepath from filekey: `local(filepath(...))`."""
    return local(storage_filepath(filekey))


@dataclass
class Settings:
    """Settings written during setup."""

    storage_dir: Union[CloudPath, Path] = None
    """Storage root. Either local dir, ``s3://bucket_name`` or ``gs://bucket_name``."""
    cache_dir: Union[Path, None] = None
    """Cache root, a local directory to cache cloud files."""
    user_name: str = None  # type: ignore
    """User name. Consider using the GitHub username."""
    user_id: Union[str, None] = None
    """User name. Consider using the GitHub username."""
    instance_name: str = None  # type: ignore
    """Name of LaminDB instance, which corresponds to exactly one backend database."""

    @property
    def cloud_storage(self) -> bool:
        """`True` if `storage_dir` is in cloud, `False` otherwise."""
        return isinstance(self.storage_dir, CloudPath)

    @property
    def _db_file(self) -> Union[Path, CloudPath]:
        """Database SQLite filepath."""
        location = self.storage_dir
        filename = str(location.stem).lower()  # type: ignore
        filepath = location / f"{filename}.lndb"  # type: ignore
        return filepath

    @property
    def db(self) -> str:
        """Database URL."""
        return f"sqlite:///{local_filepath(self._db_file)}"


def write_settings(settings: Settings):
    with open(settings_file, "w") as f:
        for key, type in get_type_hints(SettingsStore).items():
            value = getattr(settings, key)
            if value is None:
                value = "null"
            else:
                value = type(value)
            f.write(f"{key}={value}\n")


def setup_storage_dir(storage: Union[str, Path, CloudPath]) -> Union[Path, CloudPath]:
    if str(storage).startswith(("s3://", "gs://")):
        storage_dir = CloudPath(storage)
    else:
        storage_dir = Path(storage)
        if not storage_dir.exists():
            storage_dir.mkdir(parents=True)

    return storage_dir


def setup_from_store(store: SettingsStore) -> Settings:
    settings = Settings()

    settings.user_name = store.user_name
    settings.storage_dir = setup_storage_dir(store.storage_dir)
    settings.cache_dir = Path(store.cache_dir) if store.cache_dir != "null" else None
    settings.user_id = store.user_id
    settings.instance_name = store.instance_name

    return settings


def load_settings() -> Settings:
    """Return current settings."""
    if not settings_file.exists():
        logger.warning("Please setup lamindb via the CLI: lamindb setup")
        global Settings
        return Settings()
    else:
        settings_store = SettingsStore(_env_file=settings_file)
        settings = setup_from_store(settings_store)
        return settings
