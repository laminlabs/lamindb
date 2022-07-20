from pathlib import Path
from typing import Optional, Type, Union, get_type_hints  # noqa

from pydantic import BaseSettings

from ..utils.file import create_dir_if_not_exists
from .ISettings import ISettings


def create_settings_model(SettingsSchema: BaseSettings):
    class Settings(ISettings):
        def __init__(self, settings_name: str) -> None:
            self.settings_base_path = Path.home() / ".lndb"
            self.settings_name = settings_name
            self.settings_file_path = (
                self.settings_base_path / f"{self.settings_name}.env"
            )
            self.settings_schema = get_type_hints(SettingsSchema)
            self.exists = self.settings_file_path.exists()

        def setup(self, settings) -> None:
            create_dir_if_not_exists(self.settings_base_path)
            with open(self.settings_file_path, "w") as file:
                for key, value in settings.dict().items():
                    file.write(f"{key}={value}\n")
            self.exists = True

        def remove(self):
            self.settings_file_path.unlink()

        def __getitem__(self, key: str) -> Optional[str]:
            assert key in self.settings_schema, (
                f"{key} does not belong to settings {SettingsSchema} keys."
                f"Please provide one of these: {self.settings_schema.keys()}"
            )
            settings = self.__load()
            value = settings.__getattribute__(key)
            return None if value == "None" else value

        def __setitem__(self, key: str, value) -> None:
            assert key in self.settings_schema, (
                f"{key} does not belong to settings keys. "
                f"Please provide one of these: {self.settings_schema.keys()}"
            )
            assert (
                type(value) is self.settings_schema[key]
                or Optional[type(value)] is self.settings_schema[key]  # noqa
                or (  # noqa
                    isinstance(value, type(None))
                    and "Optional" in str(self.settings_schema[key])  # noqa
                )
            ), (
                f"Type {self.settings_schema[key]} expected for {key}, "
                f"but provided value {value} is {type(value)}"
            )
            settings = self.__load()
            with open(self.settings_file_path, "w") as file:
                for _key, _value in settings.dict().items():
                    if _key != key:
                        file.write(f"{_key}={_value}\n")
                file.write(f"{key}={value}\n")

        def __load(self):
            assert self.exists, (
                f"Settings file {self.settings_file_path} does not exists. "
                "Use the initmethod to create a new one."
            )
            return SettingsSchema(
                _env_file=self.settings_file_path, _env_file_encoding="utf-8"
            )

    return Settings
