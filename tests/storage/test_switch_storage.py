from pathlib import Path

import lamindb as ln
import pytest


def test_settings_switch_storage():
    ln.settings.storage = "./default_storage_unit_storage"
    assert (
        ln.settings.storage.root.resolve()
        == Path("./default_storage_unit_storage").resolve()
    )
    new_storage_location = "s3://lamindb-ci/test-settings-switch-storage"
    ln.settings.storage = new_storage_location
    assert ln.setup.settings.storage.type_is_cloud
    assert ln.setup.settings.storage.root_as_str == new_storage_location
    # root.fs contains the underlying fsspec filesystem
    # the following is set by lamindb to True for s3 by default
    assert ln.setup.settings.storage.root.fs.cache_regions
    ln.settings.storage = new_storage_location, {"cache_regions": False}
    assert not ln.setup.settings.storage.root.fs.cache_regions

    # now work with the new storage location
    new_storage = ln.Storage.get(root=new_storage_location)
    artifact = ln.Artifact(".gitignore", key="test_artifact").save()
    assert new_storage.root in artifact.path.as_posix()
    with pytest.raises(ln.setup.errors.StorageNotEmpty) as err:
        new_storage.delete()
    assert (
        "'s3://lamindb-ci/test-settings-switch-storage/.lamindb' contains 1 objects"
        in err.exconly()
    )
    artifact.delete(permanent=True)
    # switch back to default storage
    ln.settings.storage = "./default_storage_unit_storage"
