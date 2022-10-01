import shutil
from pathlib import Path
from typing import Union

import anndata
import anndata as ad
import pandas as pd
import readfcs
from anndata import AnnData
from cloudpathlib import CloudPath
from lndb_setup._settings import SettingManager

READER_FUNCS = {
    ".csv": pd.read_csv,
    ".h5ad": ad.read,
    ".feather": pd.read_feather,
    ".fcs": readfcs.read,
}


class DevFile:
    def __init__(self, settings_manager: SettingManager) -> None:
        self.settings_manager = settings_manager

    def store_file(self, localfile: Union[str, Path], storagekey: str) -> float:
        """Store arbitrary file.

        Returns size in bytes.
        """
        storagepath = self.settings_manager.instance.storage.key_to_filepath(storagekey)
        if isinstance(storagepath, CloudPath):
            storagepath.upload_from(localfile)
        else:
            shutil.copyfile(localfile, storagepath)
        size = Path(localfile).stat().st_size
        return float(size)  # because this is how we store in the db

    def delete_file(self, storagekey: str):
        """Delete arbitrary file."""
        storagepath = self.settings_manager.instance.storage.key_to_filepath(storagekey)
        storagepath.unlink()

    def load_to_memory(self, filepath: Union[str, Path]):
        filepath = Path(filepath)
        reader = READER_FUNCS.get(filepath.suffix)
        if reader is None:
            raise NotImplementedError
        else:
            return reader(filepath)

    def h5ad_to_anndata(self, filekey) -> AnnData:
        """h5ad → AnnData."""
        return anndata.read(
            self.settings_manager.instance.storage.local_filepath(filekey)
        )

    def store_png(self, filepath: str, filekey: str):
        """Store png file."""
        if not filepath.endswith("png"):
            raise ValueError()
        self.store_file(filepath, filekey)
