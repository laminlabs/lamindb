from pathlib import Path

import lamindb as ln


def test_settings_switch_storage():
    assert ln.settings.storage.resolve() == Path("./default_storage").resolve()
    ln.settings.storage = "s3://lamindb-ci"
    assert ln.setup.settings.storage.is_cloud
    assert ln.setup.settings.storage.root_as_str == "s3://lamindb-ci"
    # root.fs contains the underlying fsspec filesystem
    # the following is set by lamindb to True for s3 by default
    assert ln.setup.settings.storage.root.fs.cache_regions
    ln.settings.storage = "s3://lamindb-ci", dict(cache_regions=False)
    assert not ln.setup.settings.storage.root.fs.cache_regions
    assert ln.Storage.filter(root="s3://lamindb-ci").one_or_none() is not None
    # switch back to default storage
    ln.settings.storage = "./default_storage"
