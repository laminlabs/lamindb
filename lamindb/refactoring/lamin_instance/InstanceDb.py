from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath

from .. import model
from ..db.SQLiteLocalDbClient import SQLiteLocalDbClient


class InstanceDb:
    def __init__(self, db_client: SQLiteLocalDbClient) -> None:
        self.db_client = db_client

    # Insert user

    def insert_user_if_not_exists(self, user_email: str):
        user = self.__get_user_if_exists(user_email)
        if user is None:
            user = self.db_client.insert(model.user, {"email": user_email})
        return user

    def __get_user_if_exists(self, user_email: str):
        user = self.db_client.load_rows(model.user, {model.user.email: user_email})
        assert len(user) <= 1
        return user[0] if len(user) == 1 else None

    # Insert dobject

    def insert(
        self,
        user_id: str,
        instance_id: str,
        interface_id: str,
        dobject_id: str,
        commit_id: str,
        interface_name: str,
        path: Union[Path, CloudPath],
    ):
        interface = self.__insert_interface_if_not_exists(interface_id, interface_name)
        commit = self.__insert_commit(
            user_id,
            commit_id,
            instance_id,
            interface.id,
        )
        dobject = self.__insert_dobject(
            dobject_id,
            commit.id,
            path,
        )
        return dobject

    # Insert dobject helpers

    def __insert_interface_if_not_exists(self, interface_id: str, interface_name: str):
        interface_metadata = self.__get_interface_metadata(interface_id, interface_name)
        interface = self.__get_interface_if_exists(interface_metadata["id"])
        if interface is None:
            interface = self.db_client.insert(model.interface, interface_metadata)
        return interface

    def __insert_dobject(
        self,
        dobject_id: str,
        commit_id: str,
        path: Union[Path, CloudPath],
    ):
        dobject = self.db_client.insert(
            model.dobject,
            {
                "id": dobject_id,
                "commit_id": commit_id,
                "name": path.stem,
                "suffix": path.suffix,
            },
        )
        return dobject

    def __insert_commit(
        self, user_id: str, commit_id: str, instance_id: str, interface_id: str
    ):
        commit = self.db_client.insert(
            model.commit,
            {
                "id": commit_id,
                "user_id": user_id,
                "instance_id": instance_id,
                "interface_id": interface_id,
            },
        )
        return commit

    # insert interface helpers

    def __get_interface_if_exists(self, interface_id: str):
        interface = self.db_client.load_rows(
            model.interface, {model.interface.id: interface_id}
        )
        assert len(interface) <= 1
        return interface[0] if len(interface) == 1 else None

    def __get_interface_metadata(
        self, interface_id: str = None, interface_name: str = None
    ):
        if interface_id:
            assert interface_name is not None
            interface_metadata = self.__get_interface_other_metadata(
                interface_id, interface_name
            )
        else:
            interface_metadata = self.__get_interface_nbproject_metadata()
            assert interface_metadata["id"] is not None
            assert interface_metadata["name"] is not None
        return interface_metadata

    def __get_interface_other_metadata(self, interface_id: str, interface_name: str):
        return {
            "id": interface_id,
            "name": interface_name,
            "type": "other",
        }

    def __get_interface_nbproject_metadata(self):
        from nbproject import meta

        return {"id": meta.store.id, "name": meta.live.title, "type": "nbproject"}
