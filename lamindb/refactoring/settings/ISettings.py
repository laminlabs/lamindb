# from typing import Optional

from pydantic import BaseSettings


class ISettings:
    exists: bool

    def __init__(self, settings_name: str) -> None:
        pass

    def setup(self, settings: BaseSettings) -> None:
        pass

    # def __getitem__(self, key: str) -> Optional[str]:
    #     """Access settings by key."""
    #     pass

    # def __setitem__(self, key: str, value) -> None:
    #     """Update settings by key."""
    #     pass
