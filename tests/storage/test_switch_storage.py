from pathlib import Path

import lamindb as ln


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
    assert ln.Storage.filter(root=new_storage_location).one_or_none() is not None
    # switch back to default storage
    ln.settings.storage = "./default_storage_unit_storage"
