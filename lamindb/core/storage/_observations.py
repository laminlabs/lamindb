from __future__ import annotations

from typing import TYPE_CHECKING

import h5py
from anndata._io.specs.registry import get_spec
from lamindb_setup.core.upath import LocalPathClasses, UPath, create_mapper

from ._tiledbsoma import _open_tiledbsoma


def _X_n_obs(X):
    if "shape" in X.attrs:
        return X.attrs["shape"][0]
    else:
        return X.shape[0]


def n_observations(storepath: UPath) -> int | None:
    if storepath.is_file():
        if storepath.suffix != ".h5ad":
            return None
        with storepath.open(mode="rb") as open_obj:
            with h5py.File(open_obj, mode="r") as storage:
                return _X_n_obs(storage["X"])
    else:
        zarr_meta = {".zarray", ".zgroup"}
        tdbsoma_meta = {"__tiledb_group.tdb", "__group", "__meta"}
        is_zarr = False
        is_tdbsoma = False
        for path in storepath.iterdir():
            path_name = path.name
            if path_name in zarr_meta:
                is_zarr = True
                break
            elif path_name in tdbsoma_meta:
                is_tdbsoma = True
                break
        if is_zarr:
            try:
                import zarr
            except ImportError:
                return None
            storepath_str = storepath.as_posix()
            if isinstance(storepath, LocalPathClasses):
                open_obj = storepath_str
            else:
                open_obj = create_mapper(storepath.fs, storepath_str, check=True)
            storage = zarr.open(open_obj, mode="r")
            if get_spec(storage).encoding_type != "anndata":
                return None
            return _X_n_obs(storage["X"])
        elif is_tdbsoma:
            try:
                with _open_tiledbsoma(storepath, mode="r") as storage:
                    if "obs" in storage.keys():
                        return len(storage["obs"])
            except ImportError:
                return None
    return None
