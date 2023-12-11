from pathlib import Path

import pytest

import lamindb as ln
from lamindb.dev.storage.file import read_adata_h5ad


# https://stackoverflow.com/questions/22627659/run-code-before-and-after-each-test-in-py-test
# switch to cloud storage and back
@pytest.fixture
def switch_storage():
    cloud_storage = "s3://lamindb-ci"

    ln.settings.storage = cloud_storage

    yield cloud_storage

    ln.settings.storage = "./default_storage"


def test_local_cache():
    # check that we have local storage
    local_storage = Path("./default_storage").resolve().as_posix()
    assert ln.setup.settings.storage.root_as_str == local_storage

    test_file = ln.dev.datasets.anndata_file_pbmc68k_test()
    adata = read_adata_h5ad(test_file)

    artifact = ln.Artifact(adata, key="test_cache.h5ad")
    temp_path = artifact._local_filepath.resolve()
    assert temp_path.exists()
    assert ln.setup.settings.storage.cache_dir in temp_path.parents

    artifact.save()
    assert artifact.path.exists()
    assert not temp_path.exists()

    artifact.delete(permanent=True, storage=True)


def test_cloud_cache(switch_storage):
    # check that we have cloud storage
    assert ln.setup.settings.storage.root_as_str == switch_storage

    cache_dir = ln.setup.settings.storage.cache_dir
    assert cache_dir is not None

    test_file = ln.dev.datasets.anndata_file_pbmc68k_test()

    # test cache for saving an in-memory object
    adata = read_adata_h5ad(test_file)

    artifact = ln.Artifact(adata, key="test_cache.h5ad")
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
    artifact = ln.Artifact(test_file, key="test_cache.h5ad")
    artifact.save()
    cloud_path = artifact.path
    cache_path = ln.setup.settings.storage.cloud_to_local_no_update(cloud_path)
    assert cache_path.exists()
    assert test_file.stat().st_mtime < cache_path.stat().st_mtime
    assert cloud_path.modified.timestamp() < cache_path.stat().st_mtime

    artifact.delete(permanent=True, storage=True)
