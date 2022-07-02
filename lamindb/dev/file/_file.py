from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath, S3Client

from lamindb import setup


def filepath(filekey: Union[Path, CloudPath, str]) -> Union[Path, CloudPath]:
    """Cloud or local filepath from filekey."""
    settings = setup.settings()
    if settings.cloud_storage:
        client = S3Client(local_cache_dir=settings.cache_dir)
        return client.CloudPath(settings.storage_dir / filekey)
    else:
        return settings.storage_dir / filekey


def local(filepath: Union[Path, CloudPath]) -> Path:
    """Local (cache) filepath from filepath."""
    if setup.settings().cloud_storage:
        filepath = filepath.fspath  # type: ignore  # mypy misses CloudPath
    return filepath


def local_filepath(filekey: Union[Path, CloudPath, str]) -> Path:
    """Local (cache) filepath from filekey: `local(filepath(...))`."""
    return local(filepath(filekey))
