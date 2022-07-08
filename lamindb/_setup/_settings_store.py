from pydantic import BaseSettings


class Connector(BaseSettings):
    url: str
    key: str


class SettingsStore(BaseSettings):
    storage_dir: str
    cache_dir: str
    user_email: str
    user_id: str
    secret: str
    instance_name: str

    class Config:
        env_file = ".env"
