import shutil
from pathlib import Path

import lamindb as ln
import pytest
from lamindb.core.loaders import load_h5ad
from lamindb_setup._set_managed_storage import set_managed_storage


# https://stackoverflow.com/questions/22627659/run-code-before-and-after-each-test-in-py-test
# switch to cloud storage and back
@pytest.fixture
def switch_storage():
    cloud_storage = "s3://lamindb-ci/lamindb-unit-tests-cloud"
    set_managed_storage(cloud_storage)
    yield cloud_storage
    set_managed_storage("./default_storage_unit_core")


def test_local_cache():
    # check that we have local storage
    local_storage = Path("./default_storage_unit_core").resolve().as_posix()
    assert ln.setup.settings.storage.root_as_str == local_storage

    test_file = ln.core.datasets.anndata_file_pbmc68k_test()
    adata = load_h5ad(test_file)

    artifact = ln.Artifact.from_anndata(adata, key="test_cache.h5ad")
    temp_path = artifact._local_filepath.resolve()
    assert temp_path.exists()
    assert ln.setup.settings.storage.cache_dir in temp_path.parents

    artifact.save()
    assert artifact.path.exists()
    assert not temp_path.exists()

    artifact.delete(permanent=True, storage=True)

    # check directories
    adata_zarr_pth = Path("test_adata.zarr")
    adata.write_zarr(adata_zarr_pth)
    assert adata_zarr_pth.exists()

    artifact = ln.Artifact(adata_zarr_pth, key="test_cache.zarr")
    artifact.save()
    assert adata_zarr_pth.exists()
    assert artifact.path.exists()

    shutil.rmtree(adata_zarr_pth)
    artifact.delete(permanent=True, storage=True)

    # check directories in cache
    cache_dir = ln.setup.settings.storage.cache_dir
    adata_zarr_pth = cache_dir / "test_adata.zarr"
    adata.write_zarr(adata_zarr_pth)

    artifact = ln.Artifact(adata_zarr_pth, key="test_cache.zarr")
    assert adata_zarr_pth.exists()
    artifact.save()

    assert not adata_zarr_pth.exists()
    assert artifact.path.exists()

    artifact.delete(permanent=True, storage=True)


def test_cloud_cache(switch_storage):
    # check that we have cloud storage
    assert ln.setup.settings.storage.root_as_str == switch_storage

    cache_dir = ln.setup.settings.storage.cache_dir
    assert cache_dir is not None

    test_file = ln.core.datasets.anndata_file_pbmc68k_test()

    # test cache for saving an in-memory object
    adata = load_h5ad(test_file)

    artifact = ln.Artifact.from_anndata(adata, key="test_cache.h5ad")
    temp_path = artifact._local_filepath.resolve()
    assert cache_dir in temp_path.parents
    artifact.save()
    assert not temp_path.exists()
    cloud_path = artifact.path
    cache_path = ln.setup.settings.storage.cloud_to_local_no_update(cloud_path)
    assert cache_path.exists()
    assert cloud_path.modified.timestamp() < cache_path.stat().st_mtime

    artifact.delete(permanent=True, storage=True)

    # test cache for saving an on-disk object
    artifact = ln.Artifact.from_anndata(test_file, key="test_cache.h5ad")
    artifact.save()
    cloud_path = artifact.path
    cache_path = ln.setup.settings.storage.cloud_to_local_no_update(cloud_path)
    assert cache_path.exists()
    assert test_file.stat().st_mtime < cache_path.stat().st_mtime
    assert cloud_path.modified.timestamp() < cache_path.stat().st_mtime

    artifact.delete(permanent=True, storage=True)

    # test cache for a directory on-disk object outside the cache dir
    adata_zarr_pth = Path("test_adata.zarr")
    adata.write_zarr(adata_zarr_pth)
    artifact = ln.Artifact(adata_zarr_pth, key="test_cache.zarr")
    artifact.save()
    assert adata_zarr_pth.is_dir()
    cache_path = ln.setup.settings.storage.cloud_to_local_no_update(artifact.path)
    assert cache_path.is_dir()

    shutil.rmtree(adata_zarr_pth)
    artifact.delete(permanent=True, storage=True)

    # inside the cache dir
    adata_zarr_pth = cache_dir / "test_adata.zarr"
    adata.write_zarr(adata_zarr_pth)
    artifact = ln.Artifact(adata_zarr_pth, key="test_cache.zarr")
    assert adata_zarr_pth.exists()
    artifact.save()
    assert not adata_zarr_pth.exists()
    cache_path = ln.setup.settings.storage.cloud_to_local_no_update(artifact.path)
    assert cache_path.is_dir()

    artifact.delete(permanent=True, storage=True)
