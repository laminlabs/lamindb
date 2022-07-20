from pathlib import Path

from pydantic import BaseSettings

# user_config_dir in appdirs is weird on MacOS!
# hence, let's take home/.lndb
settings_dir = Path.home() / ".lndb"
settings_dir.mkdir(parents=True, exist_ok=True)


class instance_context:
    settings_file: Path = settings_dir / "current_instance.env"


class user_context:
    settings_file: Path = settings_dir / "current_user.env"


class InstanceSettingsStore(BaseSettings):
    storage_dir: str
    dbconfig: str

    class Config:
        env_file = ".env"


class UserSettingsStore(BaseSettings):
    user_email: str
    user_secret: str
    user_id: str

    class Config:
        env_file = ".env"


class Connector(BaseSettings):
    url: str
    key: str
