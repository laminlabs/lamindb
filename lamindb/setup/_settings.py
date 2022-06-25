import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath

root_dir = Path(__file__).parent.resolve()
settings_file = root_dir / "settings.pkl"


class description:
    storage_root = (
        "Storage root, if not a local directory, it needs to be of form"
        " `s3://bucket_name` or `gs://bucket_name`"
    )
    cache_root = "Cache root, a local directory to cache cloud files"
    user_name = "User name. Consider using the GitHub username"
    user_id = "A LaminDB user ID (8 characters, base62)"


@dataclass
class Settings:
    """Settings written during setup."""

    storage_root: Union[CloudPath, Path] = None
    cache_root: Union[Path, None] = None
    user_name: str = None  # type: ignore
    user_id: Union[str, None] = None

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
