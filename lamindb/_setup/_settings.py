from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Union, get_type_hints

from appdirs import AppDirs
from cloudpathlib import CloudPath, S3Client

from ._settings_store import (
    InstanceSettingsStore,
    UserSettingsStore,
    current_instance_settings_file,
    current_user_settings_file,
    settings_dir,
)

DIRS = AppDirs("lamindb", "laminlabs")


def storage_filepath(filekey: Union[Path, CloudPath, str]) -> Union[Path, CloudPath]:
    """Cloud or local filepath from filekey."""
    settings = load_or_create_instance_settings()
    if settings.cloud_storage:
        client = S3Client(local_cache_dir=settings.cache_dir)
        return client.CloudPath(settings.storage_dir / filekey)
    else:
        return settings.storage_dir / filekey


def cloud_to_local(filepath: Union[Path, CloudPath]) -> Path:
    """Local (cache) filepath from filepath."""
    if load_or_create_instance_settings().cloud_storage:
        filepath = filepath.fspath  # type: ignore  # mypy misses CloudPath
    Path(filepath).parent.mkdir(
        parents=True, exist_ok=True
    )  # this should not happen here
    return filepath


# conversion to Path via cloud_to_local()
# would trigger download of remote file to cache if there already
# is one
# as we don't want this, as this is a pure write operation
# we manually construct the local file path
# using the `.parts` attribute in the following line
def cloud_to_local_no_update(filepath: Union[Path, CloudPath]) -> Path:
    settings = load_or_create_instance_settings()
    if settings.cloud_storage:
        return settings.cache_dir.joinpath(*filepath.parts[1:])  # type: ignore
    return filepath


def local_filepath(filekey: Union[Path, CloudPath, str]) -> Path:
    """Local (cache) filepath from filekey: `local(filepath(...))`."""
    return cloud_to_local(storage_filepath(filekey))


# A mere tool for quick access to the docstrings above
# I thought I had it work to read from the docstrings above, but doesn't seem so
class description:
    user_email = """User email."""
    user_secret = """User login secret. Auto-generated."""
    user_id = """User ID. Auto-generated."""
    storage_dir = """Storage root. Either local dir, ``s3://bucket_name`` or ``gs://bucket_name``."""  # noqa
    _dbconfig = """Either "sqlite" or "instance_name, postgres_url"."""


def instance_from_storage(storage):
    return str(storage.stem).lower()


@dataclass
class InstanceSettings:
    """Instance settings written during setup."""

    storage_dir: Union[CloudPath, Path] = None
    """Storage root. Either local dir, ``s3://bucket_name`` or ``gs://bucket_name``."""
    _dbconfig: str = "sqlite"
    """Either "sqlite" or "instance_name, postgres_url"."""

    @property
    def cloud_storage(self) -> bool:
        """`True` if `storage_dir` is in cloud, `False` otherwise."""
        return isinstance(self.storage_dir, CloudPath)

    @property
    def cache_dir(
        self,
    ) -> Union[Path, None]:
        """Cache root, a local directory to cache cloud files."""
        if self.cloud_storage:
            cache_dir = Path(DIRS.user_cache_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)
        else:
            cache_dir = None
        return cache_dir

    @property
    def _sqlite_file(self) -> Union[Path, CloudPath]:
        """Database SQLite filepath.

        Is a CloudPath if on S3, otherwise a Path.
        """
        filename = instance_from_storage(self.storage_dir)  # type: ignore
        return storage_filepath(f"{filename}.lndb")

    @property
    def _sqlite_file_local(self):
        """If on cloud storage, update remote file."""
        return cloud_to_local_no_update(self._sqlite_file)

    def _update_cloud_sqlite_file(self):
        """If on cloud storage, update remote file."""
        if self.cloud_storage:
            sqlite_file = self._sqlite_file
            cache_file = cloud_to_local_no_update(sqlite_file)
            sqlite_file.upload_from(cache_file)

    @property
    def instance_name(self) -> str:
        """Name of LaminDB instance, which corresponds to exactly one database."""
        if self._dbconfig == "sqlite":
            return instance_from_storage(self.storage_dir)
        else:
            return self._dbconfig.split(",")[0]

    @property
    def db(self) -> str:
        """Database URL."""
        # the great thing about cloudpathlib is that it downloads the
        # remote file to cache as soon as the time stamp is out of date
        return f"sqlite:///{cloud_to_local(self._sqlite_file)}"


@dataclass
class UserSettings:
    """User Settings written during setup."""

    user_email: str = None  # type: ignore
    """User email."""
    user_secret: Union[str, None] = None
    """User login secret. Auto-generated."""
    user_id: Union[str, None] = None
    """User ID. Auto-generated."""


def write_instance_settings(settings: InstanceSettings):
    type_hints = get_type_hints(InstanceSettingsStore)
    write_settings(settings, current_instance_settings_file, type_hints)
    write_settings(settings, settings_dir / f"{settings.instance_name}.env", type_hints)


def write_user_settings(settings: UserSettings):
    type_hints = get_type_hints(UserSettingsStore)
    write_settings(settings, current_user_settings_file, type_hints)
    write_settings(settings, settings_dir / f"{settings.user_id}.env", type_hints)


def write_settings(
    settings: Union[InstanceSettings, UserSettings],
    settings_file: Path,
    type_hints: Dict[str, Any],
):
    with open(settings_file, "w") as f:
        for key, type in type_hints.items():
            settings_key = f"_{key}" if key == "dbconfig" else key
            value = getattr(settings, settings_key)
            if value is None:
                value = "null"
            else:
                value = type(value)
            f.write(f"{key}={value}\n")


def setup_storage_dir(storage: Union[str, Path, CloudPath]) -> Union[Path, CloudPath]:
    if str(storage).startswith(("s3://", "gs://")):
        storage_dir = CloudPath(storage)
    elif str(storage) == "null":
        return None
    else:
        storage_dir = Path(storage)
        storage_dir.mkdir(parents=True, exist_ok=True)
    return storage_dir


def setup_instance_from_store(store: InstanceSettingsStore) -> InstanceSettings:
    settings = InstanceSettings()
    settings.storage_dir = setup_storage_dir(store.storage_dir)
    settings._dbconfig = store.dbconfig
    return settings


def setup_user_from_store(store: UserSettingsStore) -> UserSettings:
    settings = UserSettings()
    settings.user_email = store.user_email
    settings.user_secret = store.user_secret if store.user_secret != "null" else None
    settings.user_id = store.user_id if store.user_id != "null" else None
    return settings


def load_or_create_instance_settings():
    """Return current user settings."""
    if not current_instance_settings_file.exists():
        global InstanceSettings
        return InstanceSettings()
    else:
        settings = load_instance_settings(current_instance_settings_file)
        return settings


def load_instance_settings(instance_settings_file: Path):
    settings_store = InstanceSettingsStore(_env_file=instance_settings_file)
    settings = setup_instance_from_store(settings_store)
    return settings


def load_or_create_user_settings():
    """Return current user settings."""
    if not current_user_settings_file.exists():
        global UserSettings
        return UserSettings()
    else:
        settings = load_user_settings(current_user_settings_file)
        return settings


def load_user_settings(user_settings_file: Path):
    settings_store = UserSettingsStore(_env_file=user_settings_file)
    settings = setup_user_from_store(settings_store)
    return settings
