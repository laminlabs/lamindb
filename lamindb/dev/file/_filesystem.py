from pathlib import Path
from typing import Union

import fsspec
from lndb.dev import UPath


def _infer_filesystem(path: Union[Path, UPath, str]):
    path_str = str(path)

    if isinstance(path, UPath):
        fs = path.fs
    else:
        protocol = fsspec.utils.get_protocol(path_str)
        if protocol == "s3":
            fs_kwargs = {"cache_regions": True}
        else:
            fs_kwargs = {}
        fs = fsspec.filesystem(protocol, **fs_kwargs)

    return fs, path_str
