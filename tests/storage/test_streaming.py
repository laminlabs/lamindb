import gzip
import shutil
from pathlib import Path

import anndata as ad
import h5py
import lamindb as ln
import numpy as np
import pandas as pd
import pytest
import zarr
from lamindb.core.loaders import load_h5ad
from lamindb.core.storage._anndata_accessor import _anndata_n_observations, _to_index
from lamindb.core.storage._backed_access import (
    _flat_suffixes,
    backed_access,
)
from lamindb.core.storage._polars_lazy_df import _polars_storage_options
from lamindb.core.storage._pyarrow_dataset import _open_pyarrow_dataset
from lamindb.core.storage._zarr import load_zarr
from lamindb.core.storage.objects import infer_suffix, write_to_disk


@pytest.fixture
def bad_adata_path():
    fp = ln.examples.datasets.anndata_file_pbmc68k_test()
    adata = load_h5ad(fp)
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
    test_file = ln.examples.datasets.anndata_file_pbmc68k_test()

    adata = load_h5ad(test_file)

    zarr_path = test_file.with_suffix(".zarr")
    adata.write_zarr(zarr_path)

    adata = load_zarr(zarr_path, "anndata")

    assert adata.shape == (30, 200)

    shutil.rmtree(zarr_path)


@pytest.mark.parametrize("adata_format", ["h5ad", "zarr"])
def test_backed_access(adata_format):
    fp = ln.examples.datasets.anndata_file_pbmc68k_test()
    if adata_format == "zarr":
        adata = load_h5ad(fp)

        fp = fp.with_suffix(".zarr")
        adata.write_zarr(fp)
        del adata
        # remove encoding information to check correctness of backed accessor
        store = zarr.open(fp)
        del store["obsp"]["test"].attrs["encoding-type"]
        del store["obsp"]["test"].attrs["encoding-version"]
        del store["obsm"]["X_pca"].attrs["encoding-type"]
        del store["obsm"]["X_pca"].attrs["encoding-version"]
        del store

    with pytest.raises(ValueError):
        access = backed_access(fp.with_suffix(".invalid_suffix"), using_key=None)

    # can't open anndata in write mode
    with pytest.raises(ValueError):
        access = backed_access(fp, mode="a", using_key=None)

    access = backed_access(fp, using_key=None)
    assert not access.closed

    assert isinstance(access.obs_names, pd.Index)
    assert isinstance(access.var_names, pd.Index)
    assert access.raw.shape == (30, 100)
    assert access.obsp["test"].to_memory().sum() == 30
    assert access.varp["test"].to_memory().sum() == 200
    assert access.layers["test"][0].sum() == 200

    mask = np.full(access.shape[0], False, dtype=bool)
    mask[:5] = True
    assert access[mask].X.shape == (5, 200)

    sub = access[:10]
    assert sub[:5].shape == (5, 200)
    assert sub.layers["test"].shape == sub.shape
    assert sub.raw.shape == (10, 100)
    assert sub.obsp["test"].sum() == 10
    assert sub.varp["test"].sum() == 200
    assert sub.obsm["X_pca"].shape == (10, 50)

    with pytest.raises(AttributeError):
        sub.raw.raw  # noqa: B018

    assert access[:, [1, 2, 5]].varp["test"].sum() == 3

    obs_sub = ["TCAATCACCCTTCG-8", "CGTTATACAGTACC-8", "TGCCAAGATTGTGG-7"]
    sub = access[obs_sub]
    assert sub.obs_names.tolist() == obs_sub
    assert sub.to_memory().shape == (3, 200)

    # check with a bool mask
    obs_mask = np.isin(access.obs_names, obs_sub)
    sub = access[obs_mask]
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

    with backed_access(fp, using_key=None) as access:
        assert not access.closed
        sub = access[:10]
        assert sub[:5].shape == (5, 200)
        assert sub.layers["test"].shape == sub.shape
    assert access.closed

    with backed_access(fp, using_key=None) as access:
        idx = np.array([3, 1, 2])
        assert access[:, idx].to_memory().shape == (30, 3)
        assert access[idx].to_memory().shape == (3, 200)

    if adata_format == "zarr":
        assert fp.suffix == ".zarr"
        shutil.rmtree(fp)


def test_add_column():
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test/storage"

    adata = load_h5ad(ln.examples.datasets.anndata_file_pbmc68k_test())
    zarr_path = "adata_write_mode.zarr"
    adata.write_zarr(zarr_path)

    artifact = ln.Artifact(zarr_path, description="test add_column").save()

    access = artifact.open(mode="r+")
    n_obs, n_var = access.shape
    access.add_column("obs", "ones_obs", np.ones(n_obs))
    access.add_column("var", "ones_var", np.ones(n_var))
    assert np.all(access.obs["ones_obs"] == 1)
    assert np.all(access.var["ones_var"] == 1)
    access.close()
    assert artifact.uid.endswith("0001")

    cat_col = pd.Categorical(["one"] + ["two"] * (n_obs - 1))
    with artifact.open(mode="r+") as access:
        access.add_column("obs", "cat_col", cat_col)
        assert access.obs["cat_col"].cat.categories.to_list() == ["one", "two"]
    assert artifact.uid.endswith("0002")
    # can't add in read mode
    with pytest.raises(ValueError):
        artifact.open().add_column("obs", "new_col", cat_col)

    artifact.delete(permanent=True)
    shutil.rmtree(zarr_path)

    ln.settings.storage = previous_storage


def test_to_index():
    elem_int = np.arange(3, dtype=int)
    elem_float = elem_int.astype(float)
    elem_str = elem_int.astype(str)

    assert _to_index(elem_int).dtype == "object"
    assert _to_index(elem_float).dtype == "object"
    assert _to_index(elem_str).dtype == "object"


def test_infer_suffix():
    adata = ad.AnnData()
    assert infer_suffix(adata, format="h5ad") == ".h5ad"
    with pytest.raises(ValueError):
        infer_suffix(adata, format="my format")
    with pytest.raises(NotImplementedError):
        infer_suffix(ln.Artifact)


def test_write_to_disk():
    with pytest.raises(NotImplementedError):
        write_to_disk(ln.Artifact, "path")


def test_backed_bad_format(bad_adata_path):
    access = backed_access(bad_adata_path, using_key=None)

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
    bad_adata_path.unlink()


def test_backed_zarr_not_adata():
    zarr_pth = Path("./not_adata.zarr")
    store = zarr.open(zarr_pth, mode="w")
    store["test"] = np.array(["test"])

    access = backed_access(zarr_pth)

    assert type(access).__name__ == "BackedAccessor"
    assert type(access).__module__ == "lamindb.core.storage._backed_access"
    assert access.storage["test"][...] == "test"

    shutil.rmtree(zarr_pth)


def test_anndata_open_mode():
    fp = ln.examples.datasets.anndata_file_pbmc68k_test()
    artifact = ln.Artifact(fp, key="test_adata.h5ad").save()

    with artifact.open(mode="r") as access:
        assert type(access).__name__ == "AnnDataAccessor"
        assert type(access).__module__ == "lamindb.core.storage._anndata_accessor"
    # can't open in write mode if not tiledbsoma
    with pytest.raises(ValueError):
        artifact.open(mode="w")

    artifact.delete(permanent=True, storage=True)


def test_from_lazy():
    # a different suffix in key
    with pytest.raises(ValueError):
        ln.Artifact.from_lazy(
            suffix=".zarr", overwrite_versions=True, key="mydata.h5ad"
        )

    lazy = ln.Artifact.from_lazy(
        suffix=".zarr", overwrite_versions=True, key="mydata.zarr"
    )

    store = zarr.open(lazy.path, mode="w")
    store["test"] = np.array(["test"])

    artifact = lazy.save()

    path_str = artifact.path.as_posix()
    assert ".lamindb" in path_str
    assert artifact.uid[:16] in path_str

    access = artifact.open()
    assert access.storage["test"][...] == "test"

    artifact.delete(permanent=True, storage=True)


def test_from_lazy_cloud():
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test/storage"
    lazy = ln.Artifact.from_lazy(
        suffix=".zarr", overwrite_versions=True, key="stream_test.zarr"
    )
    store = zarr.storage.FsspecStore.from_url(lazy.path.as_posix())
    group = zarr.open(store, mode="w")
    group["ones"] = np.ones(3)
    artifact = lazy.save()
    access = artifact.open()
    np.testing.assert_array_equal(access.storage["ones"][...], np.ones(3))
    artifact.delete(permanent=True, storage=True)
    ln.settings.storage = previous_storage


def test_polars_storage_options():
    storepath = ln.UPath(
        "s3://bucket/key?endpoint_url=http://localhost:9000/s3", anon=True
    )
    storage_options = _polars_storage_options(storepath)
    assert storage_options["aws_endpoint_url"] == "http://localhost:9000/s3"
    assert not storage_options["aws_virtual_hosted_style_request"]
    assert storage_options["aws_allow_http"]
    assert storage_options["aws_skip_signature"]


def test_open_dataframe_artifact():
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test/storage"
    # open from managed bucket
    artifact_remote = ln.Artifact.connect("laminlabs/lamin-dev").get(
        "iw9RRhFApeJVHC1L0001"
    )
    with artifact_remote.open(engine="polars") as ldf:
        assert ldf.collect().shape == (3, 5)

    df = pd.DataFrame({"feat1": [0, 0, 1, 1], "feat2": [6, 7, 8, 9]})
    # check as non-partitioned file
    df.to_parquet("save_df.parquet", engine="pyarrow")
    artifact_file = ln.Artifact(
        "save_df.parquet", description="Test non-partitioned parquet"
    )
    artifact_file.save()
    # cached after saving
    ds = artifact_file.open()
    assert ds.to_table().to_pandas().equals(df)
    # remove cache
    artifact_file.cache().unlink()
    # pyarrow
    ds = artifact_file.open(engine="pyarrow")
    assert ds.to_table().to_pandas().equals(df)
    # polars
    with artifact_file.open(engine="polars") as ldf:
        assert ldf.collect().to_pandas().equals(df)
    # wrong engine
    with pytest.raises(ValueError) as err:
        artifact_file.open(engine="some-other-engine")
    assert err.exconly().startswith("ValueError: Unknown engine")
    # check as partitioned folder
    df.to_parquet("save_df", engine="pyarrow", partition_cols=["feat1"])
    assert Path("save_df").is_dir()
    artifact_folder = ln.Artifact("save_df", description="Test partitioned parquet")
    artifact_folder.save()
    # cached after saving
    ds = artifact_folder.open()
    assert ds.to_table().to_pandas().equals(df[["feat2"]])
    # remove cache
    shutil.rmtree(artifact_folder.cache())
    # pyarrow
    ds = artifact_folder.open()
    assert ds.to_table().to_pandas().equals(df[["feat2"]])
    # polars
    with artifact_folder.open(engine="polars") as ldf:
        assert ldf.collect().to_pandas().equals(df[["feat2"]])
    with artifact_folder.open(engine="polars", use_fsspec=True) as ldf:
        assert ldf.collect().to_pandas().equals(df[["feat2"]])

    artifact_file.delete(permanent=True)
    artifact_folder.delete(permanent=True)

    ln.settings.storage = previous_storage


def test_open_dataframe_collection():
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test/storage"

    df = pd.DataFrame({"feat1": [0, 0, 1, 1], "feat2": [6, 7, 8, 9]})
    shard1 = ln.UPath("df1.parquet")
    shard2 = ln.UPath("df2.parquet")
    df[:2].to_parquet(shard1, engine="pyarrow")
    df[2:].to_parquet(shard2, engine="pyarrow")
    # test checking and opening local paths
    assert _flat_suffixes([shard1, ln.UPath("some.csv")]) == {".parquet", ".csv"}
    assert _open_pyarrow_dataset([shard1, shard2]).to_table().to_pandas().equals(df)

    ln.examples.datasets.file_mini_csv()

    artifact1 = ln.Artifact(shard1, key="df1.parquet").save()
    artifact2 = ln.Artifact(shard2, key="df2.parquet").save()
    artifact3 = ln.Artifact("mini.csv", key="mini.csv").save()
    artifact4 = ln.Artifact(
        "https://lamindb-test.s3.amazonaws.com/schmidt22-crispra-gws-IFNG.csv"
    ).save()

    collection1 = ln.Collection([artifact1, artifact2], key="parquet_col")
    # before saving
    # engine="pyarrow" by default
    assert collection1.open().to_table().to_pandas().equals(df)
    # after saving
    collection1.save()
    # pyarrow
    assert collection1.open(engine="pyarrow").to_table().to_pandas().equals(df)
    # polars
    with collection1.open(engine="polars") as ldf:
        assert ldf.collect().to_pandas().equals(df)
    with collection1.open(engine="polars", use_fsspec=True) as ldf:
        assert ldf.collect().to_pandas().equals(df)
    # wrong engine
    with pytest.raises(ValueError) as err:
        collection1.open(engine="some-other-engine")
    assert err.exconly().startswith("ValueError: Unknown engine")
    # different file formats
    collection2 = ln.Collection([artifact1, artifact3], key="parquet_csv_col").save()
    with pytest.raises(ValueError) as err:
        collection2.open()
    assert err.exconly().startswith(
        "ValueError: The artifacts in the collection have different file formats"
    )
    # different filesystems with pyarrow
    collection3 = ln.Collection([artifact3, artifact4], key="s3_http_col").save()
    with pytest.raises(ValueError) as err:
        collection3.open()
    assert err.exconly().startswith(
        "ValueError: The collection has artifacts with different filesystems, this is not supported"
    )

    shard1.unlink()
    shard2.unlink()

    collection1.delete(permanent=True)
    collection2.delete(permanent=True)
    collection3.delete(permanent=True)

    artifact1.delete(permanent=True)
    artifact2.delete(permanent=True)
    artifact3.delete(permanent=True)
    artifact4.delete(permanent=True, storage=False)

    Path("mini.csv").unlink(missing_ok=True)

    ln.settings.storage = previous_storage


def test_backed_wrong_suffix():
    fp = Path("test_file.txt")
    fp.write_text("test open with wrong suffix")

    artifact = ln.Artifact(fp, description="Test open wrong suffix")
    # do not save here, it just tries to open the local path
    with pytest.raises(ValueError):
        artifact.open()

    fp.unlink()


def test_anndata_n_observations(bad_adata_path):
    assert _anndata_n_observations(bad_adata_path) == 30

    assert _anndata_n_observations("./path_does_not_exist.h5ad") is None
    assert _anndata_n_observations("./path_does_not_exist.zarr") is None

    corrupted_path = Path("./corrupted.h5ad")
    shutil.copy(bad_adata_path, corrupted_path)
    with h5py.File(corrupted_path, mode="r+") as f:
        del f["obs"]
        assert "obs" not in f
    assert _anndata_n_observations(corrupted_path) is None
    corrupted_path.unlink()

    adata = ln.examples.datasets.anndata_pbmc68k_reduced()
    assert _anndata_n_observations(adata) == adata.n_obs
    zarr_path = "./test_adata_n_obs.zarr"
    adata.write_zarr(zarr_path)
    assert _anndata_n_observations(zarr_path) == adata.n_obs

    del zarr.open(zarr_path, mode="r+")["obs"].attrs["_index"]
    assert _anndata_n_observations(zarr_path) == adata.n_obs

    shutil.rmtree(zarr_path)


def _compress(input_filepath, output_filepath):
    with open(input_filepath, "rb") as f_in:
        with gzip.open(output_filepath, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


@pytest.mark.parametrize("gz_suffix", [".gz", ".tar.gz"])
def test_compressed(gz_suffix):
    adata_f = ln.examples.datasets.anndata_file_pbmc68k_test()
    adata_gz = adata_f.with_suffix(adata_f.suffix + gz_suffix)
    _compress(adata_f, adata_gz)

    artifact = ln.Artifact.from_anndata(adata_gz, key="adata.h5ad" + gz_suffix).save()
    assert artifact.n_observations == 30

    with artifact.open() as store:
        assert type(store).__name__ == "AnnDataAccessor"
        assert type(store).__module__ == "lamindb.core.storage._anndata_accessor"

    assert isinstance(artifact.load(), ad.AnnData)

    with pytest.raises(OSError):
        artifact.open(compression=None)

    artifact.delete(permanent=True)
    adata_gz.unlink()
