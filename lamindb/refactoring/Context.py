from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath

from lamindb.refactoring.instance import load_instance
from lamindb.refactoring.settings import ContextSettingsStore, create_settings_model
from lamindb.refactoring.utils.file import create_dir_if_not_exists

settings_base_path = Path(
    "/Users/fredericenard/Sources/lamin/lamindb/lamindb_data/settings"
)


class Context:
    def __init__(self, settings_base_path: Union[Path, CloudPath]) -> None:
        self.settings_base_path = settings_base_path

    def get_current_user(self):
        context_settings = self.__create_or_get_context_settings()
        assert context_settings["current_user_email"] is not None
        assert context_settings["current_user_secret"] is not None
        return {
            "current_user_email": context_settings["current_user_email"],
            "current_user_secret": context_settings["current_user_secret"],
        }

    def set_current_user(self, email: str, secret: str):
        context_settings = self.__create_or_get_context_settings()
        context_settings["current_user_email"] = email
        context_settings["current_user_secret"] = secret

    def get_current_instance(self):
        context_settings = self.__create_or_get_context_settings()
        assert context_settings["current_instance_name"] is not None
        assert context_settings["current_user_email"] is not None
        assert context_settings["current_user_secret"] is not None
        current_instance = load_instance(
            context_settings["current_instance_name"],
            self.settings_base_path,
            context_settings["current_user_email"],
        )
        return current_instance

    def set_current_instance(self, instance_name: str):
        context_settings = self.__create_or_get_context_settings()
        context_settings["current_instance_name"] = instance_name

    def __create_or_get_context_settings(self):
        create_dir_if_not_exists(self.settings_base_path)
        ContextSettings = create_settings_model(ContextSettingsStore)
        settings = ContextSettings("context", settings_base_path)
        settings_store = ContextSettingsStore()
        if not settings.exists:
            settings.setup(settings_store)
        return settings


context = Context(settings_base_path)
