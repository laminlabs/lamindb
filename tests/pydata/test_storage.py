import concurrent.futures
from uuid import UUID

import lamindb as ln
import pytest
from lamindb.models import SQLRecord


# we need this test both in the core and the storage/cloud tests
# because the internal logic that retrieves information about other instances
# depends on whether the current instance is managed on the hub
def test_reference_storage_location(ccaplog):
    ln.Artifact("s3://lamindata/iris_studies/study0_raw_images")
    assert ln.Storage.get(root="s3://lamindata").instance_uid == "4XIuR0tvaiXM"
    assert "storage location s3://lamindata is already marked" in ccaplog.text


def test_create_storage_locations_parallel():
    root: str = "nonregistered_storage"

    def create_storage() -> str:
        ln.Storage(root=root).save()  # type: ignore
        return root

    n_parallel = 3
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_parallel) as executor:
        futures = [executor.submit(create_storage) for i in range(n_parallel)]
        _ = [future.result() for future in concurrent.futures.as_completed(futures)]

    storage = ln.Storage.get(root__endswith=root)
    storage.delete()


def test_storage_host_property(tmp_path):
    local = ln.Storage(root=(tmp_path / "host-local").as_posix(), host="test-host")
    assert local.host == "test-host"

    cloud = ln.Storage(root="s3://lamindb-ci/test-host-property", type="s3")
    assert cloud.host is None


def test_save_deletes_hub_record_when_local_save_fails(tmp_path, monkeypatch):
    current_instance_uid = ln.setup.settings.instance.uid

    class _Settings:
        root_as_str = (tmp_path / "rollback-hub").as_posix()
        instance_uid = current_instance_uid
        type = "local"
        region = None
        _uid = "abcdefghij12"
        _uuid = UUID(int=1)

    deleted = []
    monkeypatch.setattr(
        "lamindb.models.storage.init_storage",
        lambda *args, **kwargs: (_Settings(), "hub-record-created"),
    )
    monkeypatch.setattr(
        "lamindb.models.storage.delete_storage_record",
        lambda ssettings: deleted.append(ssettings),
    )

    def fail_save(self, *args, **kwargs):
        raise RuntimeError("local save failed")

    monkeypatch.setattr(SQLRecord, "save", fail_save)

    storage = ln.Storage(root=_Settings.root_as_str, instance_uid=current_instance_uid)
    with pytest.raises(RuntimeError, match="local save failed"):
        storage.save()
    assert len(deleted) == 1
    assert deleted[0]._uuid == _Settings._uuid
    assert storage._created_hub_record is False
