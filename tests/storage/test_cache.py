import shutil
from pathlib import Path
from time import sleep

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
    set_managed_storage("./default_storage_unit_storage")


def test_local_cache():
    # check that we have local storage
    local_storage = Path("./default_storage_unit_storage").resolve().as_posix()
    assert ln.setup.settings.storage.root_as_str == local_storage

    test_file = ln.core.datasets.anndata_file_pbmc68k_test()
    adata = load_h5ad(test_file)

    artifact = ln.Artifact.from_anndata(adata, key="test_cache.h5ad")
    temp_path = artifact._local_filepath.resolve()
    assert temp_path.exists()
    assert ln.setup.settings.cache_dir in temp_path.parents

    artifact.save()
    assert artifact.path.exists()
    assert not temp_path.exists()

    artifact.delete(permanent=True)

    # check directories
    adata_zarr_pth = Path("test_adata.zarr")
    adata.write_zarr(adata_zarr_pth)
    assert adata_zarr_pth.exists()

    artifact = ln.Artifact(adata_zarr_pth, key="test_cache.zarr").save()
    assert adata_zarr_pth.exists()
    assert artifact.path.exists()
    assert artifact.path.name != artifact.key

    shutil.rmtree(adata_zarr_pth)
    artifact.delete(permanent=True)

    # check directories in cache
    cache_dir = ln.setup.settings.cache_dir
    adata_zarr_pth = cache_dir / "test_adata.zarr"
    adata.write_zarr(adata_zarr_pth)

    artifact = ln.Artifact(adata_zarr_pth, key="test_cache.zarr")
    assert adata_zarr_pth.exists()
    artifact.save()

    assert not adata_zarr_pth.exists()
    assert artifact.path.exists()
    assert artifact.path.name != artifact.key

    artifact.delete(permanent=True)


def test_cloud_cache(switch_storage):
    # check that we have cloud storage
    assert ln.setup.settings.storage.root_as_str == switch_storage

    cache_dir = ln.setup.settings.cache_dir
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
    cache_path = artifact._cache_path
    assert cache_path.exists()
    assert (
        cache_path == cache_dir / "lamindb-ci/lamindb-unit-tests-cloud/test_cache.h5ad"
    )
    assert cloud_path.modified.timestamp() < cache_path.stat().st_mtime

    artifact.delete(permanent=True)

    # test cache for saving an on-disk object
    artifact = ln.Artifact.from_anndata(test_file, key="test_cache.h5ad")
    artifact.save()
    cloud_path = artifact.path
    cache_path = artifact._cache_path
    assert cache_path.exists()
    assert (
        cache_path == cache_dir / "lamindb-ci/lamindb-unit-tests-cloud/test_cache.h5ad"
    )
    assert test_file.stat().st_mtime < cache_path.stat().st_mtime
    assert cloud_path.modified.timestamp() < cache_path.stat().st_mtime

    artifact.delete(permanent=True)

    # test cache for a directory on-disk object outside the cache dir
    adata_zarr_pth = Path("test_adata.zarr")
    adata.write_zarr(adata_zarr_pth)
    artifact = ln.Artifact(adata_zarr_pth, key="test_cache.zarr")
    artifact.save()
    assert adata_zarr_pth.is_dir()
    cache_path = artifact._cache_path
    assert cache_path.is_dir()
    assert (
        cache_path == cache_dir / "lamindb-ci/lamindb-unit-tests-cloud/test_cache.zarr"
    )

    shutil.rmtree(adata_zarr_pth)
    artifact.delete(permanent=True)

    # inside the cache dir
    adata_zarr_pth = cache_dir / "test_adata.zarr"
    adata.write_zarr(adata_zarr_pth)
    artifact = ln.Artifact(adata_zarr_pth, key="test_cache.zarr")
    assert adata_zarr_pth.exists()
    artifact.save()
    assert not adata_zarr_pth.exists()
    cache_path = artifact._cache_path
    assert cache_path.is_dir()
    assert (
        cache_path == cache_dir / "lamindb-ci/lamindb-unit-tests-cloud/test_cache.zarr"
    )

    artifact.delete(permanent=True)


def test_cloud_cache_versions(switch_storage):
    adata = load_h5ad(ln.core.datasets.anndata_file_pbmc68k_test())

    cache_dir = ln.setup.settings.cache_dir
    assert cache_dir is not None

    artifact = ln.Artifact.from_anndata(adata, key="test_cache.h5ad")
    assert ln.settings.cache_dir in artifact._local_filepath.parents
    artifact.save()
    cache_path_v1 = artifact.cache()
    assert cache_path_v1.exists()
    assert (
        cache_path_v1
        == cache_dir / "lamindb-ci/lamindb-unit-tests-cloud/test_cache.h5ad"
    )
    cache_path_v1.unlink()
    artifact.cache(print_progress=False)
    assert cache_path_v1.exists()
    assert (
        cache_path_v1
        == cache_dir / "lamindb-ci/lamindb-unit-tests-cloud/test_cache.h5ad"
    )
    timestamp_v1 = cache_path_v1.stat().st_mtime
    # hope it is enough to avoid random timestamp problems further
    sleep(1)
    # new version
    adata.obs["test_cache"] = "test"
    artifact_v2 = ln.Artifact.from_anndata(
        adata, key="test_cache.h5ad", revises=artifact
    )
    assert ln.settings.cache_dir in artifact_v2._local_filepath.parents
    artifact_v2.save()
    assert artifact_v2.is_latest
    assert not artifact.is_latest
    cache_path_v2 = artifact_v2.cache()
    assert cache_path_v2.exists()
    assert (
        cache_path_v2
        == cache_dir / "lamindb-ci/lamindb-unit-tests-cloud/test_cache.h5ad"
    )
    assert cache_path_v2.stat().st_mtime > timestamp_v1
    cache_path_v2.unlink()
    artifact_v2.cache(mute=True)
    assert cache_path_v2.exists()
    assert (
        cache_path_v2
        == cache_dir / "lamindb-ci/lamindb-unit-tests-cloud/test_cache.h5ad"
    )
    assert "test_cache" in load_h5ad(cache_path_v2).obs.columns
    cache_mtime = cache_path_v2.stat().st_mtime
    assert cache_mtime == artifact_v2.path.modified.timestamp()
    assert cache_mtime > timestamp_v1
    # old version cache ignores key
    cache_path_v1 = artifact.cache()
    assert cache_path_v1.exists()
    assert cache_path_v1.name == f"{artifact.uid}.h5ad"

    artifact_v2.versions.delete(permanent=True)


def test_corrupted_cache_local():
    filepath = ln.core.datasets.anndata_file_pbmc68k_test()
    artifact = ln.Artifact.from_anndata(filepath, key="test_corrupt_cache_local.h5ad")
    artifact.save()
    # corrupt cache
    with open(artifact._cache_path, "r+b") as f:
        f.write(b"corruption")
    # just raises an exception, nothing to re-sync on local
    with pytest.raises(OSError):
        artifact.load()
    with pytest.raises(OSError):
        artifact.open()

    artifact.delete(permanent=True)


def test_corrupted_cache_cloud(switch_storage):
    # check that we have cloud storage
    assert ln.setup.settings.storage.root_as_str == switch_storage

    filepath = ln.core.datasets.anndata_file_pbmc68k_test()
    artifact = ln.Artifact.from_anndata(filepath, key="test_corrupt_cache_cloud.h5ad")
    artifact.save()
    # corrupt cache
    # sleep not to reset cache mtime to a smaller value
    # it is increased artificially on cache copying in save
    # so due to lower granularity of cloud mtimes and fast code execution
    # after the change cache mtime can become smaller than cloud mtime
    sleep(1)
    with open(artifact._cache_path, "r+b") as f:
        f.write(b"corruption")
    assert artifact._cache_path.stat().st_mtime > artifact.path.stat().st_mtime
    # check that it is indeed corrupted
    with pytest.raises(OSError):
        load_h5ad(artifact.cache())
    # should load successfully
    artifact.load()
    # check open also
    assert artifact._cache_path.exists()
    with open(artifact._cache_path, "r+b") as f:
        f.write(b"corruption")
    # should open successfully
    with artifact.open():
        pass
    # corrupted cache has been deleted
    assert not artifact._cache_path.exists()

    artifact.delete(permanent=True)
