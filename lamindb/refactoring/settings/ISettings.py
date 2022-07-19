from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath
from pydantic import BaseSettings
from pyparsing import Optional


class ISettings:
    exists: bool

    def __init__(
        self, settings_name: str, settings_base_path: Union[Path, CloudPath]
    ) -> None:
        pass

    def setup(self, settings: BaseSettings) -> None:
        pass

    def __getitem__(self, key: str) -> Optional[str]:
        """Access settings by key."""
        pass

    def __setitem__(self, key: str, value) -> None:
        """Update settings by key."""
        pass
