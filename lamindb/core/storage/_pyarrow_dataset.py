from __future__ import annotations

from typing import TYPE_CHECKING

import pyarrow.dataset
from lamindb_setup.core.upath import LocalPathClasses

if TYPE_CHECKING:
    from pyarrow.dataset import Dataset as PyArrowDataset
    from upath import UPath


PYARROW_SUFFIXES = (".parquet", ".csv", ".json", ".orc", ".arrow", ".feather", ".ipc")


def _is_pyarrow_dataset(paths: UPath | list[UPath]) -> bool:
    # it is assumed here that the paths exist
    # we don't check here that the filesystem is the same
    # but this is a requirement for pyarrow.dataset.dataset
    if isinstance(paths, list):
        path_list = paths
    elif paths.is_file():
        path_list = [paths]
    else:
        path_list = [path for path in paths.rglob("*") if path.suffix != ""]
    suffix = None
    for path in path_list:
        # ignore .gz, only check the real suffix
        path_suffixes = path.suffixes
        path_suffix = (
            path_suffixes[-2]
            if len(path_suffixes) > 1 and ".gz" in path_suffixes
            else path.suffix
        )
        if path_suffix not in PYARROW_SUFFIXES:
            return False
        elif suffix is None:
            suffix = path_suffix
        elif path_suffix != suffix:
            return False
    return True


def _open_pyarrow_dataset(paths: UPath | list[UPath]) -> PyArrowDataset:
    if isinstance(paths, list):
        path0 = paths[0]
        if isinstance(path0, LocalPathClasses):
            paths_str, filesystem = [path.as_posix() for path in paths], None
        else:
            paths_str, filesystem = [path.path for path in paths], path0.fs
    elif isinstance(paths, LocalPathClasses):
        paths_str, filesystem = paths.as_posix(), None
    else:
        paths_str, filesystem = paths.path, paths.fs

    return pyarrow.dataset.dataset(paths_str, filesystem=filesystem)
