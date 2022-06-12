from pathlib import Path
from typing import Union

import anndata
from anndata._core.anndata import AnnData
from cloudpathlib import CloudPath, S3Client
from typeguard import typechecked

from ...dev import logger


class File:
    """Access a file.

    Params
    ------
    path
        File path.
    """

    def __init__(self, path: Union[Path, str]) -> None:
        # global variables
        from lamindb._configuration import cache_root  # isort: skip
        from lamindb._configuration import cloud_storage, storage_root

        self._cloud_storage = cloud_storage
        self._storage_root = (
            storage_root  # not a Path! would strip double slash in s3://!
        )
        self._cache_root = Path(cache_root)
        # the path
        if cloud_storage:
            client = S3Client(local_cache_dir=cache_root)
            self._path = client.CloudPath(
                f"{self._storage_root}/{path}"
            )  # cannot use Path here!
        else:
            self._path = Path(self._storage_root) / path

    @property
    def path(self) -> Union[Path, CloudPath]:
        """The file path."""
        return self._path


class h5ad(File):
    """Save or load an `.h5ad` file.

    Params
    ------
    path
        File path.
    """

    def read(self) -> AnnData:
        """Load file to object."""
        path = self.path
        if self._cloud_storage:
            path = self.path.fspath  # type: ignore  # mypy misses CloudPath
        return anndata.read(path)

    @typechecked
    def write(self, adata: AnnData) -> None:
        """Save object.

        Params
        ------
        adata
           AnnData object.
        """
        if self._cloud_storage:
            # conversion to Path would trigger download of cache file below
            # hence, we just use the `.name` attribute in the following line
            cache_file = self._cache_root.joinpath(*self.path.parts[1:])
            if not cache_file.parent.exists():
                cache_file.parent.mkdir()
            logger.debug(f"writing cache file: {cache_file}")
            adata.write(cache_file)
            logger.debug("uploading cache file")
            self.path.upload_from(cache_file)  # type: ignore  # mypy misses CloudPath
            # in principle, we could write the cache file to disk again so that
            # the time stamp is newer than the one in the cloud, avoiding
            # download to access the just written cache however, cloudpath lib
            # complains about the newer cache file and will attempt download,
            # currently there doesn't seem to be a solution for this
            # adata.write(cache_file)
        else:
            adata.write(self.path)
