from pathlib import Path

from lamindb.dev.storage import delete_storage
from lamindb.dev.storage._file import read_adata_h5ad
from lamindb.dev.storage._zarr import read_adata_zarr, write_adata_zarr


def test_anndata_io():
    adata = read_adata_h5ad("s3://lamindb-ci/lndb-storage/pbmc68k.h5ad")

    def callback(*args, **kwargs):
        pass

    zarr_path = Path("./pbmc68k.zarr")
    write_adata_zarr(adata, zarr_path, callback)

    adata = read_adata_zarr(zarr_path)

    assert adata.shape == (70, 765)

    delete_storage(zarr_path)
