from pathlib import Path

from anndata import AnnData
from lamin_logger import logger
from lndb_setup import settings
from typeguard import typechecked


@typechecked
def anndata_to_h5ad(adata: AnnData, filekey: str) -> Path:
    """AnnData → h5ad."""
    path = settings.instance.storage.key_to_filepath(filekey)
    if settings.instance.cloud_storage:
        cache_file = settings.instance.storage.cloud_to_local_no_update(path)  # type: ignore  # noqa
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
