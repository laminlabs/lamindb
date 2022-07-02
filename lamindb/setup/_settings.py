from dataclasses import dataclass
from pathlib import Path
from typing import Union, get_type_hints

from cloudpathlib import CloudPath
from pydantic import BaseSettings

root_dir = Path(__file__).parent.resolve()
settings_file = root_dir / ".env"


@dataclass
class Settings:
    """Settings written during setup."""

    storage_root: Union[CloudPath, Path] = None
    """Storage root. Either local dir, ``s3://bucket_name`` or ``gs://bucket_name``."""
    cache_dir: Union[Path, None] = None
    """Cache root, a local directory to cache cloud files."""
    user_name: str = None  # type: ignore
    """User name. Consider using the GitHub username."""
    user_id: Union[str, None] = None
    """User name. Consider using the GitHub username."""

    @property
    def cloud_storage(self) -> bool:
        """`True` if `storage_root` is in cloud, `False` otherwise."""
        return isinstance(self.storage_root, CloudPath)

    @property
    def _db_file(self) -> Path:
        """Database SQLite filepath."""
        if not self.cloud_storage:
            location = self.storage_root
        else:
            location = self.cache_dir
        filename = str(location.stem).lower()  # type: ignore
        filepath = location / f"{filename}.lndb"  # type: ignore
        return filepath

    @property
    def db(self) -> str:
        """Database URL."""
        return f"sqlite:///{self._db_file}"


# A mere tool for quick access to the docstrings above
# I thought I had it work to read from the docstrings above, but doesn't seem so
class description:
    storage_root = """Storage root. Either local dir, ``s3://bucket_name`` or ``gs://bucket_name``."""  # noqa
    cache_dir = """Cache root, a local directory to cache cloud files."""
    user_name = """User name. Consider using the GitHub username."""
    user_id = """User name. Consider using the GitHub username."""


class SettingsStore(BaseSettings):
    storage_root: str
    cache_dir: str
    user_name: str
    user_id: str

    class Config:
        env_file = ".env"


def _write(settings: Settings):
    with open(settings_file, "w") as f:
        for key, type in get_type_hints(SettingsStore).items():
            value = getattr(settings, key)
            if value is None:
                value = "null"
            else:
                value = type(value)
            f.write(f"{key}={value}\n")
