from __future__ import annotations

from typing import TYPE_CHECKING

import pyarrow.dataset
from lamindb_setup.core.upath import LocalPathClasses

if TYPE_CHECKING:
    from upath import UPath


PYARROW_SUFFIXES = (".parquet", ".csv", ".json", ".orc", ".arrow", ".feather")


def _is_pyarrow_dataset(paths: UPath | list[UPath]) -> bool:
    # it is assumed here that the paths exist
    if isinstance(paths, list):
        suffixes = {path.suffix for path in paths}
    elif paths.is_file():
        suffixes = {paths.suffix}
    else:
        suffixes = {path.suffix for path in paths.rglob("*") if path.suffix != ""}
    return len(suffixes) == 1 and suffixes.pop() in PYARROW_SUFFIXES


def _open_pyarrow_dataset(paths: UPath | list[UPath]) -> pyarrow.dataset.Dataset:
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
