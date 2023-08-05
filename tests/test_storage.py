from pathlib import Path

import pandas as pd
import pytest

from lamindb.dev.storage import delete_storage
from lamindb.dev.storage._backed_access import backed_access
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


@pytest.mark.parametrize("adata_format", ["h5ad", "zarr"])
def test_backed_access(adata_format):
    fp = Path("tests/test-files/pbmc68k.h5ad")
    if adata_format == "zarr":
        adata = read_adata_h5ad(fp)

        def callback(*args, **kwargs):
            pass

        fp = fp.with_suffix(".zarr")
        write_adata_zarr(adata, fp, callback)
        del adata

    with pytest.raises(ValueError):
        access = backed_access(fp.with_suffix(".invalid_suffix"))

    access = backed_access(fp)

    assert isinstance(access.obs_names, pd.Index)
    assert isinstance(access.var_names, pd.Index)
    assert access.raw.shape == (30, 100)
    assert access.obsp["test"].to_memory().sum() == 30
    assert access.varp["test"].to_memory().sum() == 200
    assert access.layers["test"][0].sum() == 200

    sub = access[:10]
    assert sub[:5].shape == (5, 200)
    assert sub.layers["test"].shape == sub.shape
    assert sub.raw.shape == (10, 100)
    assert sub.obsp["test"].sum() == 10
    assert sub.varp["test"].sum() == 200

    with pytest.raises(AttributeError):
        sub.raw.raw

    assert access[:, [1, 2, 5]].varp["test"].sum() == 3

    obs_sub = ["TCAATCACCCTTCG-8", "CGTTATACAGTACC-8", "TGCCAAGATTGTGG-7"]
    sub = access[obs_sub]
    assert sub.obs_names.tolist() == obs_sub

    var_sub = ["SSU72", "PARK7", "RBP7"]
    sub = access[:, var_sub]
    assert sub.var_names.tolist() == var_sub

    assert access.to_memory().shape == (30, 200)
    assert sub.to_memory().shape == (30, 3)

    if adata_format == "zarr":
        assert fp.suffix == ".zarr"
        delete_storage(fp)
