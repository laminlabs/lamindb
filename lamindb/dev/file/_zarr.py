import fsspec
from anndata import AnnData
from anndata._io import read_zarr, write_zarr


def read_adata_zarr(filepath) -> AnnData:
    store = fsspec.get_mapper(filepath, check=True)
    adata = read_zarr(store)

    return adata


def write_adata_zarr(adata, filepath, **kwargs):
    # todo: delete if exists already
    store = fsspec.get_mapper(filepath, create=True)
    write_zarr(store, adata, **kwargs)
