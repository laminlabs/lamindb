import os
from pathlib import Path

from anndata import AnnData
from lamin_logger import logger
from lndb import settings
from typeguard import typechecked


@typechecked
def anndata_to_h5ad(adata: AnnData, filekey: str) -> Path:
    """AnnData → h5ad."""
    path = settings.instance.storage.key_to_filepath(filekey)
    if settings.instance.storage.is_cloud:
        cache_file = settings.instance.storage.cloud_to_local_no_update(path)  # type: ignore  # noqa
        cache_file.parent.mkdir(exist_ok=True)
        logger.debug(f"Writing cache file: {cache_file}.")
        adata.write(cache_file)
        logger.debug("Uploading cache file.")
        path.upload_from(cache_file)  # type: ignore
        # to avoid download from the cloud within synchronization
        mtime = path.modified.timestamp()
        os.utime(cache_file, times=(mtime, mtime))
    else:
        adata.write(path)
        cache_file = path
    return cache_file
