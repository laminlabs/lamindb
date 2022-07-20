from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath

from ..instance.CommitManager import CommitManager
from ..instance.InstanceDb import InstanceDb


class Do:
    def __init__(
        self,
        instance_id: str,
        instance_db: InstanceDb,
        commit_manager: CommitManager,
    ) -> None:
        self.instance_id = instance_id
        self.instance_db = instance_db
        self.commit_manager = commit_manager

    # Load

    def load(self, entity):
        return self.instance_db.db_client.load_table(entity)

    # Add

    def add(
        self,
        path: Union[str, Path, CloudPath],
        interface_id: str = None,
        interface_name: str = None,
    ):
        path = self.__cast_path_if_str(path)
        self.commit_manager.add(
            path,
            interface_id,
            interface_name,
        )

    # Commit

    def commit(self):
        self.commit_manager.commit()

    # Status

    def status(self):
        return self.commit_manager.added

    # Utils

    def __cast_path_if_str(self, path: Union[str, Path, CloudPath]):
        if type(path) is str:
            return Path(path)
        else:
            return path
