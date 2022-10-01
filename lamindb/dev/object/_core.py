from pathlib import Path

from anndata import AnnData
from lamin_logger import logger
from lndb_setup._settings import SettingManager
from pandas import DataFrame
from typeguard import typechecked


class DevObject:
    def __init__(self, settings_manager: SettingManager) -> None:
        self.settings_manager = settings_manager

    def infer_suffix(self, dmem):
        """Infer LaminDB storage file suffix from a data object."""
        if isinstance(dmem, AnnData):
            return ".h5ad"
        elif isinstance(dmem, DataFrame):
            return ".feather"
        else:
            raise NotImplementedError

    def write_to_file(self, dmem, filepath: str):
        if isinstance(dmem, AnnData):
            dmem.write(filepath)
        elif isinstance(dmem, DataFrame):
            try:
                dmem.to_feather(filepath)
            except ValueError:
                dmem.reset_index().to_feather(filepath)
        else:
            raise NotImplementedError

    @typechecked
    def anndata_to_h5ad(self, adata: AnnData, filekey: str) -> Path:
        """AnnData → h5ad."""
        path = self.settings_manager.instance.storage.key_to_filepath(filekey)
        if self.settings_manager.instance.cloud_storage:
            cache_file = self.settings_manager.instance.storage.cloud_to_local_no_update(path)  # type: ignore  # noqa
            cache_file.parent.mkdir(exist_ok=True)
            logger.debug(f"Writing cache file: {cache_file}.")
            adata.write(cache_file)
            logger.debug("Uploading cache file.")
            path.upload_from(cache_file)  # type: ignore  # mypy misses CloudPath
            # In principle, we could write the cache file to disk again so that
            # the time stamp is newer than the one in the cloud, avoiding
            # download to access the just written cache. However, cloudpathlib
            # complains about the newer cache file and will attempt download,
            # currently there doesn't seem to be a solution for this
        else:
            adata.write(path)
            cache_file = path
        return cache_file
