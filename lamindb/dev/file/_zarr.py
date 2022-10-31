import warnings
from typing import Optional

import fsspec
import scipy.sparse as sparse
import zarr
from anndata import AnnData
from anndata._io import read_zarr
from anndata._io.specs import write_elem

from ..object._anndata_sizes import _size_elem, _size_raw, size_adata


def read_adata_zarr(storepath) -> AnnData:
    if not isinstance(storepath, str):
        storepath = str(storepath)
    store = fsspec.get_mapper(storepath, check=True)
    adata = read_zarr(store)

    return adata


def write_adata_zarr(
    adata: AnnData, storepath, callback=None, chunks=None, **dataset_kwargs
):
    if not isinstance(storepath, str):
        storepath = str(storepath)

    store = fsspec.get_mapper(storepath, create=True)
    f = zarr.open(store, mode="w")

    adata.strings_to_categoricals()
    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)

    f.attrs.setdefault("encoding-type", "anndata")
    f.attrs.setdefault("encoding-version", "0.1.0")

    adata_size = None
    cumulative_val = 0

    def _cb(key_write: Optional[str] = None):
        nonlocal adata_size
        nonlocal cumulative_val

        if callback is None:
            return None
        if adata_size is None:
            adata_size = size_adata(adata)
        if key_write is None:
            # begin or finish
            if cumulative_val < adata_size:
                callback(adata_size, adata_size if cumulative_val > 0 else 0)
            return None

        elem = getattr(adata, key_write, None)
        if elem is None:
            return None
        elem_size = _size_raw(elem) if key_write == "raw" else _size_elem(elem)
        if elem_size == 0:
            return None

        cumulative_val += elem_size
        callback(adata_size, cumulative_val)

    def _write_elem_cb(f, k, elem, dataset_kwargs):
        write_elem(f, k, elem, dataset_kwargs=dataset_kwargs)
        _cb(k)

    _cb(None)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="zarr")

        if chunks is not None and not isinstance(adata.X, sparse.spmatrix):
            _write_elem_cb(
                f, "X", adata.X, dataset_kwargs=dict(chunks=chunks, **dataset_kwargs)
            )
        else:
            _write_elem_cb(f, "X", adata.X, dataset_kwargs=dataset_kwargs)
        for elem in ("obs", "var"):
            _write_elem_cb(f, elem, getattr(adata, elem), dataset_kwargs=dataset_kwargs)
        for elem in ("obsm", "varm", "obsp", "varp", "layers", "uns"):
            _write_elem_cb(
                f, elem, dict(getattr(adata, elem)), dataset_kwargs=dataset_kwargs
            )
        _write_elem_cb(f, "raw", adata.raw, dataset_kwargs=dataset_kwargs)
    # todo: fix size less than total at the end
    _cb(None)
