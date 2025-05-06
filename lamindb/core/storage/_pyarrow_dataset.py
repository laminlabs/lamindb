from __future__ import annotations

from typing import TYPE_CHECKING

import pyarrow.dataset
from lamindb_setup.core.upath import LocalPathClasses

if TYPE_CHECKING:
    from pyarrow.dataset import Dataset as PyArrowDataset
    from upath import UPath


PYARROW_SUFFIXES = (".parquet", ".csv", ".json", ".orc", ".arrow", ".feather", ".ipc")


def _open_pyarrow_dataset(paths: UPath | list[UPath], **kwargs) -> PyArrowDataset:
    if isinstance(paths, list):
        # a single path can be a directory, but a list of paths
        # has to be a flat list of files
        paths_str = []
        path0 = paths[0]
        if isinstance(path0, LocalPathClasses):
            path_to_str = lambda p: p.as_posix()
            filesystem = None
        else:
            path_to_str = lambda p: p.path
            filesystem = path0.fs
        for path in paths:
            if (
                getattr(path, "protocol", None) not in {"http", "https"}
                and path.is_dir()
            ):
                paths_str += [path_to_str(p) for p in path.rglob("*") if p.suffix != ""]
            else:
                paths_str.append(path_to_str(path))
    elif isinstance(paths, LocalPathClasses):
        paths_str, filesystem = paths.as_posix(), None
    else:
        paths_str, filesystem = paths.path, paths.fs

    return pyarrow.dataset.dataset(paths_str, filesystem=filesystem, **kwargs)
