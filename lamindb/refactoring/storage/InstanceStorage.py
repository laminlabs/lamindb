import shutil
from pathlib import Path
from typing import Union

from anndata import AnnData
from cloudpathlib import CloudPath


class InstanceStorage:
    def __init__(self, storage_path: Union[Path, CloudPath]) -> None:
        self.storage_path = storage_path

    def save(self, to_save: Union[AnnData, Path, CloudPath], key: str):
        if isinstance(to_save, AnnData):
            self.save_anndata(to_save, key)
        else:
            self.save_file_from_path(to_save, key)

    def save_anndata(self, anndata: AnnData, key: str):
        to_path = self.storage_path / key
        anndata.write(to_path)

    def save_file_from_path(self, from_path: Union[Path, CloudPath], key: str):
        to_path = self.storage_path / key
        if isinstance(to_path, CloudPath):
            self.__save_cloud(from_path, to_path)
        else:
            self.__save_local(from_path, to_path)

    def delete(self, key: str):
        to_path = self.storage_path / key
        to_path.unlink()

    def __save_local(self, from_path: Path, to_path: Path):
        shutil.copyfile(from_path, to_path)

    def __save_cloud(self, from_path: CloudPath, to_path: CloudPath):
        to_path.upload_from(from_path)
