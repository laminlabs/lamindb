from pydantic import BaseSettings


class SettingsStore(BaseSettings):
    storage_dir: str
    cache_dir: str
    user_name: str
    user_id: str
    instance_name: str

    class Config:
        env_file = ".env"
