import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath

root_dir = Path(__file__).parent.resolve()
settings_file = root_dir / "settings.pkl"


@dataclass
class Settings:
    """Settings written during setup."""

    storage_root: Union[CloudPath, Path] = None
    """Storage root. Either local dir, ``s3://bucket_name`` or ``gs://bucket_name``."""
    cache_root: Union[Path, None] = None
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
            location = self.cache_root
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
    cache_root = """Cache root, a local directory to cache cloud files."""
    user_name = """User name. Consider using the GitHub username."""
    user_id = """User name. Consider using the GitHub username."""


def _write(settings: Settings):
    with open(settings_file, "wb") as f:
        pickle.dump(settings, f, protocol=4)


def settings() -> Settings:
    """Return current settings."""
    if not settings_file.exists():
        print("WARNING: Please setup lamindb via the CLI: lamindb setup")
        global Settings
        return Settings()
    else:
        with open(settings_file, "rb") as f:
            settings = pickle.load(f)
        return settings
