import shutil

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr

import lamindb as ln
from lamindb.dev.storage import delete_storage
from lamindb.dev.storage._backed_access import backed_access
from lamindb.dev.storage._zarr import read_adata_zarr, write_adata_zarr
from lamindb.dev.storage.file import read_adata_h5ad
from lamindb.dev.storage.object import infer_suffix, write_to_file


@pytest.fixture
def bad_adata_path():
    fp = ln.dev.datasets.anndata_file_pbmc68k_test()
    adata = read_adata_h5ad(fp)
    to = fp.with_name("pbmc68k_bad.h5ad")
    shutil.copy(fp, to)
    fp = to
    file = h5py.File(fp, mode="r+")
    for field_name in ("obs", "var"):
        field = getattr(adata, field_name).to_records()
        formats = []
        for name, (dt, _) in field.dtype.fields.items():
            if dt == "O":
                new_dt = str(field[name].astype(str).dtype).replace("<U", "S")
            else:
                new_dt = dt
            formats.append((name, new_dt))
        del file[field_name]
        file.create_dataset(field_name, data=field.astype(formats))
    del file["X"].attrs["encoding-type"]
    del file["X"].attrs["encoding-version"]
    del file["obsp"]["test"].attrs["encoding-type"]
    del file["obsp"]["test"].attrs["encoding-version"]
    file.close()
    return fp


def test_anndata_io():
    test_file = ln.dev.datasets.anndata_file_pbmc68k_test()

    adata = read_adata_h5ad(test_file)

    def callback(*args, **kwargs):
        pass

    zarr_path = test_file.with_suffix(".zarr")
    write_adata_zarr(adata, zarr_path, callback)

    adata = read_adata_zarr(zarr_path)

    assert adata.shape == (30, 200)

    delete_storage(zarr_path)


@pytest.mark.parametrize("adata_format", ["h5ad", "zarr"])
def test_backed_access(adata_format):
    fp = ln.dev.datasets.anndata_file_pbmc68k_test()
    if adata_format == "zarr":
        adata = read_adata_h5ad(fp)

        def callback(*args, **kwargs):
            pass

        fp = fp.with_suffix(".zarr")
        write_adata_zarr(adata, fp, callback)
        del adata
        # remove encoding information to check correctness of backed accessor
        store = zarr.open(fp)
        del store["obsp"]["test"].attrs["encoding-type"]
        del store["obsp"]["test"].attrs["encoding-version"]
        del store["obsm"]["X_pca"].attrs["encoding-type"]
        del store["obsm"]["X_pca"].attrs["encoding-version"]
        del store

    with pytest.raises(ValueError):
        access = backed_access(fp.with_suffix(".invalid_suffix"))

    access = backed_access(fp)
    assert not access.closed

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
    assert sub.obsm["X_pca"].shape == (10, 50)

    with pytest.raises(AttributeError):
        sub.raw.raw

    assert access[:, [1, 2, 5]].varp["test"].sum() == 3

    obs_sub = ["TCAATCACCCTTCG-8", "CGTTATACAGTACC-8", "TGCCAAGATTGTGG-7"]
    sub = access[obs_sub]
    assert sub.obs_names.tolist() == obs_sub
    assert sub.to_memory().shape == (3, 200)

    idx = np.array([1, 2, 5])
    sub = access[idx]
    assert sub.raw.shape == (3, 100)
    assert sub.to_memory().shape == (3, 200)

    var_sub = ["SSU72", "PARK7", "RBP7"]
    sub = access[:, var_sub]
    assert sub.var_names.tolist() == var_sub

    assert access.to_memory().shape == (30, 200)
    assert sub.to_memory().shape == (30, 3)

    access.close()
    assert access.closed
    del access

    with backed_access(fp) as access:
        assert not access.closed
        sub = access[:10]
        assert sub[:5].shape == (5, 200)
        assert sub.layers["test"].shape == sub.shape
    assert access.closed

    if adata_format == "zarr":
        assert fp.suffix == ".zarr"
        delete_storage(fp)


def test_infer_suffix():
    import anndata as ad

    adata = ad.AnnData()
    assert infer_suffix(adata, adata_format="h5ad") == ".h5ad"
    with pytest.raises(ValueError):
        infer_suffix(adata, adata_format="my format")
    with pytest.raises(NotImplementedError):
        infer_suffix(ln.File)


def test_write_to_file():
    with pytest.raises(NotImplementedError):
        write_to_file(ln.File, "path")


def test_backed_bad_format(bad_adata_path):
    access = backed_access(bad_adata_path)

    assert access.obsp["test"].to_memory().sum() == 30

    sub = access[:10]

    assert sub.X.shape == (10, 200)
    assert sub.obsp["test"].sum() == 10

    assert isinstance(sub.obs, pd.DataFrame)
    assert isinstance(sub.var, pd.DataFrame)
    assert isinstance(sub.obs_names, pd.Index)
    assert isinstance(sub.var_names, pd.Index)

    assert sub.to_memory().shape == (10, 200)

    access.close()
    delete_storage(bad_adata_path)
