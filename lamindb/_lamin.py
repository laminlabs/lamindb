from lamindb.dev import Dev
from lamindb.schema import Schema
from lndb_setup._settings import SettingManager


class Lamin:

    def __init__(self, settings_manager: SettingManager) -> None:
        self.settings_manager = settings_manager
        self.schema = Schema(self.settings_manager)
        self.dev = Dev(self.settings_manager)
