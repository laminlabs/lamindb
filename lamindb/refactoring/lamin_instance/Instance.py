from ..db.SQLiteLocalDbClient import SQLiteLocalDbClient
from ..settings.create_settings_model import create_settings_model
from ..settings.settings_store import InstanceSettingsStore
from ..storage.InstanceStorage import InstanceStorage
from .CommitManager import CommitManager
from .Do import Do
from .InstanceDb import InstanceDb
from .Schema import Schema


class Instance:
    def __init__(self, settings, user_email) -> None:
        self.settings = settings
        self.instance_id = self.settings["instance_id"]
        self.instance_name = self.settings["instance_name"]
        self.instance_db = self.__load_instance_db()
        self.instance_storage = self.__load_instance_storage()
        self.user = self.__create_user_if_not_exists(user_email)
        self.commit_manager = self.__create_commit_manager()
        self.do = Do(
            instance_id=self.instance_id,
            instance_db=self.instance_db,
            commit_manager=self.commit_manager,
        )
        self.schema = Schema(self.instance_db.db_client.engine)

    def __create_user_if_not_exists(self, user_email):
        return self.instance_db.insert_user_if_not_exists(user_email)

    def __create_commit_manager(self):
        commit_manager = CommitManager(
            user_id=self.user.id,
            instance_id=self.instance_id,
            instance_db=self.instance_db,
            instance_storage=self.instance_storage,
        )
        return commit_manager

    def __load_instance_db(self):
        if self.settings["db_config"] == "sqlite":
            db_client = SQLiteLocalDbClient(
                instance_name=self.settings["instance_name"],
                db_base_path=self.settings["db_base_path"],
            )
        else:
            raise NotImplementedError()
        instance_db = InstanceDb(db_client)
        return instance_db

    def __load_instance_storage(self):
        return InstanceStorage(
            self.settings["storage_base_path"] / self.settings["instance_name"]
        )


def load_instance(instance_name, user_email):
    InstanceSettings = create_settings_model(InstanceSettingsStore)
    instance_settings = InstanceSettings(instance_name)
    instance = Instance(settings=instance_settings, user_email=user_email)
    return instance
