from pathlib import Path

import lamindb as ln
import pytest
from lamindb_setup.core._hub_core import get_storage_records_for_instance


def check_storage_location_on_hub_exists(uid: str):
    all_storage_records = get_storage_records_for_instance(
        ln.setup.settings.instance._id
    )
    length = len([r for r in all_storage_records if r["lnid"] == uid])
    if length not in {0, 1}:
        raise AssertionError(
            f"Expected 0 or 1 storage records for uid {uid}, found {length}."
        )
    return length == 1


def test_reference_storage_location(ccaplog):
    ln.Artifact("s3://lamindata/iris_studies/study0_raw_images")
    assert ln.Storage.get(root="s3://lamindata").instance_uid == "4XIuR0tvaiXM"
    # assert (
    #     "referenced read-only storage location at s3://lamindata, is managed by instance with uid 4XIuR0tvaiXM"
    #     in ccaplog.text
    # )


def test_switch_delete_storage_location():
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
    assert ln.setup.settings.storage.root.exists()

    # now work with the new storage location
    new_storage = ln.Storage.get(root=new_storage_location)
    assert check_storage_location_on_hub_exists(new_storage.uid)
    artifact = ln.Artifact(".gitignore", key="test_artifact").save()
    assert new_storage.root in artifact.path.as_posix()

    # artifacts exist
    with pytest.raises(AssertionError) as err:
        new_storage.delete()
    assert "Cannot delete storage holding artifacts." in err.exconly()

    artifact.delete(permanent=True, storage=False)
    # still some files in there
    with pytest.raises(ln.setup.errors.StorageNotEmpty) as err:
        new_storage.delete()
    assert (
        "'s3://lamindb-ci/test-settings-switch-storage/.lamindb' contains 1 objects"
        in err.exconly()
    )

    # now delete the artifact so that the storage location is empty
    artifact.path.unlink()
    with pytest.raises(AssertionError) as err:
        new_storage.delete()
    assert (
        "Cannot delete the current storage location, switch to another."
        in err.exconly()
    )

    # check all attempts unsuccessful so far
    assert check_storage_location_on_hub_exists(new_storage.uid)

    # switch back to default storage
    ln.settings.storage = "./default_storage_unit_storage"
    storage_marker = ln.UPath(new_storage_location) / ".lamindb/storage_uid.txt"
    assert storage_marker.exists()
    new_storage.delete()
    assert not check_storage_location_on_hub_exists(new_storage.uid)
    assert not storage_marker.exists()
