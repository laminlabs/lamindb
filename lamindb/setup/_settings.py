from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath
from pydantic import BaseModel, Field

root_dir = Path(__file__).parent.resolve()
settings_file = root_dir / "settings.json"


class description:
    storage_root = (
        "Storage root, if not a local directory, it needs to be of form"
        " 's3://bucket_name' or 'gs://bucket_name'"
    )
    cache_root = "Cache root, a local directory to cache cloud files."
    user_name = "User name. Consider using the GitHub username."
    user_id = "A LaminDB user ID (8 characters, base62)."


class Settings(BaseModel):
    """Settings written during setup."""

    # pydantic does not understand CloudPath... hence we need this str field
    storage_root_str: str = Field(description=description.storage_root)
    """{description.storage_root}."""
    cache_root: Union[Path, None] = Field(description=description.cache_root)
    """{description.cache_root}."""
    user_name: str = Field(description=description.user_name)
    """{description.user_name}."""
    user_id: Union[str, None] = Field(default=None, description=description.user_id)
    """{description.user_id}."""

    @property
    def storage_root(self) -> Union[Path, CloudPath]:
        """{description.storage_root}."""
        if self.storage_root_str.startswith(("s3://", "gs://")):
            return CloudPath(self.storage_root_str)
        else:
            return Path(self.storage_root_str)

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

    def _write(self):
        with open(settings_file, "w") as f:
            f.write(self.json())


def load():
    with open(settings_file) as f:
        settings_json = f.read()
    return Settings.parse_raw(settings_json)
