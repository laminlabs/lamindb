import warnings
from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath

# from lamindb.refactoring.settings.ISettings import ISettings
from lamindb.refactoring.utils.file import create_dir_if_not_exists
from lamindb.refactoring.utils.id import id_instance

from .Context import Context
from .settings import InstanceSettingsStore, create_settings_model


class SetupInstance:
    @staticmethod
    def setup_if_not_exists(
        instance_name: str,
        db_base_path: Union[Path, CloudPath],
        storage_base_path: Union[Path, CloudPath],
        db_config: str = "sqlite",
    ) -> None:
        instance_id = id_instance()
        instance_settings = SetupInstance.__setup_instance_settings_if_not_exists(
            instance_id,
            instance_name,
            db_base_path,
            storage_base_path,
            db_config,
        )
        SetupInstance.__setup_instance_db_if_not_exists(instance_settings)
        SetupInstance.__setup_instance_storage_if_not_exists(instance_settings)
        Context.set_current_instance(instance_name)

    @staticmethod
    def remove_instance(instance_name: str):
        InstanceSettings = create_settings_model(InstanceSettingsStore)
        settings = InstanceSettings(instance_name)
        settings["db_base_path"].unlink()
        (settings["storage_base_path"] / instance_name).rmdir()

    @staticmethod
    def __setup_instance_settings_if_not_exists(
        instance_id: str,
        instance_name: str,
        db_base_path: Union[Path, CloudPath],
        storage_base_path: Union[Path, CloudPath],
        db_config: str = "sqlite",
    ):
        InstanceSettings = create_settings_model(InstanceSettingsStore)
        settings = InstanceSettings(instance_name)
        if not settings.exists:
            settings_store = InstanceSettingsStore(
                instance_id=instance_id,
                instance_name=instance_name,
                db_base_path=db_base_path.absolute(),
                storage_base_path=storage_base_path.absolute(),
                db_config=db_config,
            )
            settings.setup(settings_store)
        else:
            warnings.warn(f"Instance {instance_name} already installed.")
        return settings

    @staticmethod
    def __setup_instance_db_if_not_exists(settings) -> None:
        create_dir_if_not_exists(settings["db_base_path"])

    @staticmethod
    def __setup_instance_storage_if_not_exists(settings):
        create_dir_if_not_exists(
            settings["storage_base_path"] / settings["instance_name"]
        )
