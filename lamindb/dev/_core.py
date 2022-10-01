from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath
from lamindb import db
from lndb_setup._settings import SettingManager
from lnschema_core import dobject, type
from .object._core import DevObject
from .file._file import DevFile
from .db._core import DevDb


class Dev:
    def __init__(self, settings_manager: SettingManager) -> None:
        self.settings_manager = settings_manager
        self.object = DevObject(settings_manager)
        self.file = DevFile(settings_manager)
        self.db = DevDb(settings_manager)

    def get_name_suffix_from_filepath(self, filepath: Union[Path, CloudPath]):
        suffix = "".join(filepath.suffixes)
        name = filepath.name.replace(suffix, "")
        return name, suffix

    def storage_key_from_dobject(self, dobj: dobject):
        return f"{dobj.id}{dobj.suffix}"

    def storage_key_from_triple(self, dobj_id: str, dobj_suffix: str):
        return f"{dobj_id}{dobj_suffix}"

    def filepath_from_dobject(self, dobj: dobject):
        storage_key = self.storage_key_from_dobject(dobj)
        filepath = self.settings_manager.instance.storage.key_to_filepath(storage_key)
        return filepath

    def track_usage(self, dobject_id, usage_type: type.usage):
        usage_id = getattr(db.insert, "usage")(
            type=usage_type,
            user_id=self.settings_manager.user.id,
            dobject_id=dobject_id,
        )

        return usage_id
