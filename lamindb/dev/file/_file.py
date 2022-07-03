import shutil
from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath, S3Client

from lamindb import setup


def storage_filepath(filekey: Union[Path, CloudPath, str]) -> Union[Path, CloudPath]:
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
    return local(storage_filepath(filekey))


def store_file(filepath: Union[str, Path], filekey: str):
    """Store png file."""
    storage_path = storage_filepath(filekey)
    if isinstance(storage_path, CloudPath):
        storage_path.upload_from(filepath)
    else:
        shutil.copyfile(filepath, storage_path)
