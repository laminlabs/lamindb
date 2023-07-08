from pathlib import Path

import pytest

import lamindb as ln
from lamindb.dev.storage.file import read_adata_h5ad


# https://stackoverflow.com/questions/22627659/run-code-before-and-after-each-test-in-py-test
# switch to cloud storage and back
# needed to isolate this test from the others
@pytest.fixture(autouse=True)
def switch_storage():
    ln.settings.storage = "s3://lamindb-ci"

    yield

    ln.settings.storage = "./default_storage"


def test_cache():
    cache_dir = ln.setup.settings.storage.cache_dir
    assert cache_dir is not None

    test_file = Path("tests/test-files/pbmc68k.h5ad")

    # test cache for saving an in-memory object
    adata = read_adata_h5ad(test_file)

    file = ln.File(adata, key="test_cache.h5ad")
    temp_path = file._local_filepath.resolve()
    assert cache_dir in temp_path.parents
    file.save()
    assert not temp_path.exists()
    cloud_path = file.path()
    cache_path = ln.setup.settings.storage.cloud_to_local_no_update(cloud_path)
    assert cache_path.exists()
    assert cloud_path.modified.timestamp() < cache_path.stat().st_mtime

    file.delete(storage=True)

    # test cache for saving an on-disk object
    file = ln.File(test_file, key="test_cache.h5ad")
    file.save()
    cloud_path = file.path()
    cache_path = ln.setup.settings.storage.cloud_to_local_no_update(cloud_path)
    assert cache_path.exists()
    assert test_file.stat().st_mtime < cache_path.stat().st_mtime
    assert cloud_path.modified.timestamp() < cache_path.stat().st_mtime

    file.delete(storage=True)
