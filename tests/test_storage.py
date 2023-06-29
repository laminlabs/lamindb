from pathlib import Path

import h5py

from lamindb.dev.storage import _infer_filesystem, delete_storage
from lamindb.dev.storage._backed_access import AnnDataAccessor
from lamindb.dev.storage._zarr import read_adata_zarr, write_adata_zarr
from lamindb.dev.storage.file import read_adata_h5ad


def test_anndata_io():
    test_files = Path("tests/test-files")

    adata = read_adata_h5ad(test_files / "pbmc68k.h5ad")

    def callback(*args, **kwargs):
        pass

    zarr_path = test_files / "pbmc68k.zarr"
    write_adata_zarr(adata, zarr_path, callback)

    adata = read_adata_zarr(zarr_path)

    assert adata.shape == (30, 200)

    delete_storage(zarr_path)


def test_backed_access():
    fs, file_path_str = _infer_filesystem("tests/test-files/pbmc68k.h5ad")
    conn = fs.open(file_path_str, mode="rb")
    storage = h5py.File(conn, mode="r")

    access = AnnDataAccessor(conn, storage, "pbmc68k.h5ad")

    assert access.raw.shape == (30, 100)

    sub = access[:10]
    assert sub[:5].shape == (5, 200)
    assert sub.raw.shape == (10, 100)
