import hashlib
from pathlib import Path
from typing import Dict, Union

from cloudpathlib import CloudPath

from lamindb.refactoring.storage.InstanceStorage import InstanceStorage
from lamindb.refactoring.utils.id import id_commit

from .InstanceDb import InstanceDb


class CommitManager:
    def __init__(
        self,
        user_id: str,
        instance_id: str,
        instance_db: InstanceDb,
        instance_storage: InstanceStorage,
    ) -> None:
        self.user_id = user_id
        self.instance_id = instance_id
        self.instance_db = instance_db
        self.instance_storage = instance_storage
        self.commit_id = id_commit()
        self.added: Dict[str, Dict] = {}

    def add(
        self,
        path: Union[Path, CloudPath],
        interface_id: str = None,
        interface_name: str = None,
    ):
        checksum = hashlib.md5(open(path, "rb").read()).hexdigest()
        self.added[checksum] = {
            "user_id": self.user_id,
            "instance_id": self.instance_id,
            "interface_id": interface_id,
            "commit_id": self.commit_id,
            "dobject_id": checksum,
            "path": path,
            "interface_name": interface_name,
        }

    def commit(self):
        for id, parameters in self.added.items():
            self.instance_db.insert(**parameters)
            self.instance_storage.save(
                parameters["path"],
                f"""{parameters["dobject_id"]}{parameters["path"].suffix}""",
            )
        self.commit_id = id_commit()
        self.added = {}

    def status(self):
        return self.added
