import concurrent.futures
from uuid import UUID

import lamindb as ln
import pytest
from lamindb.models import SQLRecord
from lamindb_setup.core._settings_storage import mark_storage_root_file
from lamindb_setup.core.upath import UPath


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


def _hub_settings(storage_root: UPath, current_instance_uid: str):
    class _Settings:
        root = storage_root
        root_as_str = storage_root.as_posix()
        instance_uid = current_instance_uid
        type = "local"
        region = None
        _uid = "abcdefghij12"
        _uuid = UUID(int=1)

    return _Settings


def _write_root_mark(root: UPath) -> UPath:
    mark_file = mark_storage_root_file(root)
    mark_file.parent.mkdir(parents=True, exist_ok=True)
    mark_file.write_text("abcdefghij12\n")
    return mark_file


def test_save_deletes_hub_record_when_local_save_fails(tmp_path, monkeypatch):
    current_instance_uid = ln.setup.settings.instance.uid
    root = UPath(tmp_path / "rollback-hub")
    settings_cls = _hub_settings(root, current_instance_uid)
    mark_file = _write_root_mark(root)
    deleted = []
    monkeypatch.setattr(
        "lamindb.models.storage.init_storage",
        lambda *args, **kwargs: (settings_cls(), "hub-record-created"),
    )
    monkeypatch.setattr(
        "lamindb.models.storage.delete_storage_record",
        lambda ssettings: deleted.append(ssettings),
    )

    def fail_save(self, *args, **kwargs):
        raise RuntimeError("local save failed")

    monkeypatch.setattr(SQLRecord, "save", fail_save)

    storage = ln.Storage(root=root.as_posix(), instance_uid=current_instance_uid)
    with pytest.raises(RuntimeError, match="local save failed"):
        storage.save()
    assert len(deleted) == 1
    assert deleted[0]._uuid == settings_cls._uuid
    assert storage._created_hub_record is False
    assert not mark_file.exists()


def test_save_keeps_hub_record_when_local_insert_completed(tmp_path, monkeypatch):
    current_instance_uid = ln.setup.settings.instance.uid
    root = UPath(tmp_path / "keep-hub")
    settings_cls = _hub_settings(root, current_instance_uid)
    mark_file = _write_root_mark(root)
    deleted = []
    monkeypatch.setattr(
        "lamindb.models.storage.init_storage",
        lambda *args, **kwargs: (settings_cls(), "hub-record-created"),
    )
    monkeypatch.setattr(
        "lamindb.models.storage.delete_storage_record",
        lambda ssettings: deleted.append(ssettings),
    )

    def fail_after_insert(self, *args, **kwargs):
        self._state.adding = False
        raise RuntimeError("post save failed")

    monkeypatch.setattr(SQLRecord, "save", fail_after_insert)

    storage = ln.Storage(root=root.as_posix(), instance_uid=current_instance_uid)
    with pytest.raises(RuntimeError, match="post save failed"):
        storage.save()
    assert deleted == []
    assert storage._created_hub_record is True
    assert mark_file.exists()


def test_save_uses_fields_edited_before_first_save(tmp_path, monkeypatch):
    original = (tmp_path / "original-root").as_posix()
    edited_path = tmp_path / "edited-root"
    edited_path.mkdir()
    edited = edited_path.as_posix()
    seen = {}

    def fake_init_storage(root, **kwargs):
        seen["root"] = root
        seen["region"] = kwargs["region"]

        class _Settings:
            root_as_str = root
            instance_uid = ln.setup.settings.instance.uid
            type = "s3"
            region = kwargs["region"]
            _uid = "abcdefghij12"
            _uuid = None

        return _Settings(), "hub-record-not-created"

    monkeypatch.setattr("lamindb.models.storage.init_storage", fake_init_storage)
    storage = ln.Storage(root=original, type="local", region="host-a")
    storage.root = edited
    storage.type = "s3"
    storage.region = "host-b"
    storage.save()
    assert seen["root"] == edited
    assert seen["region"] == "host-b"
    assert storage.root == edited
    assert storage.type == "s3"
    assert storage.region == "host-b"
    storage.delete()


def test_deferred_storage_space_change_updates_hub(tmp_path, monkeypatch):
    root = (tmp_path / "deferred-space").as_posix()
    storage = ln.Storage(root=root).save()
    space = ln.Space(name=f"storage-space-{storage.uid}").save()
    updates = []
    monkeypatch.setattr(
        "lamindb.models.storage.update_storage_with_space",
        lambda **kwargs: updates.append(kwargs),
    )
    deferred = ln.Storage.objects.only("id", "uid").get(id=storage.id)
    deferred.space = space
    deferred.save()
    assert updates == [{"storage_lnid": storage.uid, "space_lnid": space.uid}]
    storage.delete()
    space.delete()
