from pathlib import Path
from typing import Optional, Union

from cloudpathlib import CloudPath
from pydantic import BaseSettings


class InstanceSettingsStore(BaseSettings):
    instance_id: str
    instance_name: str
    db_base_path: Union[Path, CloudPath]
    storage_base_path: Union[Path, CloudPath]
    db_config: str


class ContextSettingsStore(BaseSettings):
    current_instance_name: Optional[str]
    current_user_email: Optional[str]
    current_user_secret: Optional[str]
