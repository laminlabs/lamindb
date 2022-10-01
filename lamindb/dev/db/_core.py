import sqlmodel as sqm
from lndb_setup._settings import SettingManager


class DevDb:
    def __init__(self, settings_manager: SettingManager) -> None:
        self.settings_manager = settings_manager

    def session(self) -> sqm.Session:
        """Connection session to DB engine.

        Returns a `sqlmodel.Session` object.
        """
        return sqm.Session(self.settings_manager.instance.db_engine())
