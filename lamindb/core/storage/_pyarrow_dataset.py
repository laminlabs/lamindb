from __future__ import annotations

from typing import TYPE_CHECKING

import pyarrow.dataset
from lamindb_setup.core.upath import LocalPathClasses

if TYPE_CHECKING:
    from upath import UPath


PYARROW_SUFFIXES = (".parquet", ".csv", ".json", ".orc", ".arrow", ".feather")


def _is_pyarrow_dataset(path: UPath) -> bool:
    # it is assumed here that path exists
    if path.is_file():
        return path.suffix in PYARROW_SUFFIXES
    else:
        objects = path.rglob("*")
        suffixes = {object.suffix for object in objects if object.suffix != ""}
        return len(suffixes) == 1 and suffixes.pop() in PYARROW_SUFFIXES


def _open_pyarrow_dataset(path: UPath) -> pyarrow.dataset.Dataset:
    if isinstance(path, LocalPathClasses):
        path_str, _filesytem = path.as_posix(), None
    else:
        path_str, filesystem = path.path, path.fs

    return pyarrow.dataset.dataset(path_str, filesystem=filesystem)
