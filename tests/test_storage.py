import shutil
from pathlib import Path

import h5py
import lamindb as ln
import numpy as np
import pandas as pd
import pytest
import tiledbsoma
import tiledbsoma.io
import zarr
from lamindb.core.storage._backed_access import (
    AnnDataAccessor,
    BackedAccessor,
    backed_access,
)
from lamindb.core.storage._tiledbsoma import create_tiledbsoma_store
from lamindb.core.storage._zarr import read_adata_zarr, write_adata_zarr
from lamindb.core.storage.objects import infer_suffix, write_to_disk
from lamindb.core.storage.paths import read_adata_h5ad


@pytest.fixture
def bad_adata_path():
    fp = ln.core.datasets.anndata_file_pbmc68k_test()
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
    test_file = ln.core.datasets.anndata_file_pbmc68k_test()

    adata = read_adata_h5ad(test_file)

    def callback(*args, **kwargs):
        pass

    zarr_path = test_file.with_suffix(".zarr")
    write_adata_zarr(adata, zarr_path, callback)

    adata = read_adata_zarr(zarr_path)

    assert adata.shape == (30, 200)

    shutil.rmtree(zarr_path)


@pytest.mark.parametrize("adata_format", ["h5ad", "zarr"])
def test_backed_access(adata_format):
    fp = ln.core.datasets.anndata_file_pbmc68k_test()
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


def test_infer_suffix():
    import anndata as ad

    adata = ad.AnnData()
    assert infer_suffix(adata, adata_format="h5ad") == ".h5ad"
    with pytest.raises(ValueError):
        infer_suffix(adata, adata_format="my format")
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
    store["test"] = "test"

    access = backed_access(zarr_pth)

    assert isinstance(access, BackedAccessor)
    assert access.storage["test"][...] == "test"

    shutil.rmtree(zarr_pth)


def test_anndata_open_mode():
    fp = ln.core.datasets.anndata_file_pbmc68k_test()
    artifact = ln.Artifact(fp, key="test_adata.h5ad")
    artifact.save()

    with artifact.open(mode="r") as access:
        assert isinstance(access, AnnDataAccessor)
    # can't open in write mode if not tiledbsoma
    with pytest.raises(ValueError):
        artifact.open(mode="w")

    artifact.delete(permanent=True, storage=True)


@pytest.mark.parametrize("storage", [None, "s3://lamindb-test"])
def test_backed_tiledbsoma(storage):
    if storage is not None:
        previous_storage = ln.setup.settings.storage.root_as_str
        ln.settings.storage = "s3://lamindb-test"

    test_file = ln.core.datasets.anndata_file_pbmc68k_test()
    tiledbsoma.io.from_h5ad("test.tiledbsoma", test_file, "RNA")

    artifact_soma = ln.Artifact("test.tiledbsoma", description="test tiledbsoma")
    artifact_soma.save()

    # copied to cache on .save()
    cache_path = artifact_soma.cache()

    hash_on_disk = artifact_soma.hash
    with artifact_soma.open(mode="w") as store:
        assert store.__class__.__name__ == "ExperimentTrack"
    if storage is not None:
        # hash in the cloud will be different from hash on disk
        # therefore the artifact will be updated
        assert artifact_soma.hash != hash_on_disk
        # delete the cached store on hash change
        assert not cache_path.exists()
    else:
        # hash stays the same
        assert artifact_soma.hash == hash_on_disk
        assert artifact_soma.path == cache_path

    experiment = artifact_soma.open()  # mode="r" by default
    assert isinstance(experiment, tiledbsoma.Experiment)
    experiment.close()

    # wrong mode, should be either r or w for tiledbsoma
    with pytest.raises(ValueError):
        artifact_soma.open(mode="p")

    # run deprecated backed
    with artifact_soma.backed():
        pass

    artifact_soma.versions.delete(permanent=True, storage=True)
    shutil.rmtree("test.tiledbsoma")

    if storage is not None:
        ln.settings.storage = previous_storage


@pytest.mark.parametrize("storage", [None, "s3://lamindb-test"])
def test_create_tiledbsoma_store(storage):
    if storage is not None:
        previous_storage = ln.setup.settings.storage.root_as_str
        ln.settings.storage = "s3://lamindb-test"

    fp = ln.core.datasets.anndata_file_pbmc68k_test()

    transform = ln.Transform(name="test create tiledbsoma store")
    transform.save()
    run = ln.Run(transform)
    run.save()

    experiment_path = (ln.settings.storage.root / "test_create.tiledbsoma").as_posix()
    artifact = create_tiledbsoma_store(experiment_path, fp, run, measurement_name="RNA")
    artifact.save()

    with artifact.open() as store:
        assert "lamin_run_uid" in store["obs"].schema.names

    artifact.delete(permanent=True, storage=True)

    if storage is not None:
        ln.settings.storage = previous_storage
