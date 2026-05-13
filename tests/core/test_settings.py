import lamindb as ln
import pytest


def test_settings_repr():
    repr_str = repr(ln.settings)

    lines = repr_str.split("\n")
    assert "Settings" in lines[0]
    assert all(line.startswith("  ") for line in lines[1:])

    content = "\n".join(lines[1:])
    assert content.find("instance:") < content.find("storage:")
    assert content.find("storage:") < content.find("verbosity:")
    assert content.find("verbosity:") < content.find("track_run_inputs:")


def test_storage_setter_raises_on_foreign_managed_storage(tmp_path):
    storage = ln.Storage(root=(tmp_path / "foreign-managed-storage").as_posix()).save()
    storage.instance_uid = "_not_exists_"
    storage.save()

    with pytest.raises(ValueError) as error:
        ln.settings.storage = storage.root
    assert (
        error.exconly()
        == f"ValueError: Storage '{storage.root}' exists in another instance (_not_exists_), cannot write to it from here."
    )
    storage.delete()


def test_local_storage_setter_raises_on_foreign_managed_storage(tmp_path):
    storage = ln.Storage(
        root=(tmp_path / "foreign-managed-local-storage").as_posix()
    ).save()
    storage.instance_uid = "_not_exists_"
    storage.save()

    with pytest.raises(ValueError) as error:
        ln.settings.local_storage = storage.root
    assert (
        error.exconly()
        == f"ValueError: Storage '{storage.root}' exists in another instance (_not_exists_), cannot write to it from here."
    )
    storage.delete()


def test_storage_setter_raises_on_unmanaged_storage(tmp_path):
    storage = ln.Storage(root=(tmp_path / "unmanaged-storage").as_posix()).save()
    storage.instance_uid = None
    storage.save()

    with pytest.raises(ValueError) as error:
        ln.settings.storage = storage.root
    assert (
        error.exconly()
        == f"ValueError: Storage '{storage.root}' is not managed by any instance, cannot write to it from here."
    )
    storage.delete()


def test_local_storage_setter_raises_on_unmanaged_storage(tmp_path):
    storage = ln.Storage(root=(tmp_path / "unmanaged-local-storage").as_posix()).save()
    storage.instance_uid = None
    storage.save()

    with pytest.raises(ValueError) as error:
        ln.settings.local_storage = storage.root
    assert (
        error.exconly()
        == f"ValueError: Storage '{storage.root}' is not managed by any instance, cannot write to it from here."
    )
    storage.delete()
