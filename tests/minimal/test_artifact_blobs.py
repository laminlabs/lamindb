"""Constructor tests for blob artifacts.

These run in the minimal suite, which does not install numpy. Constructing an
unstructured artifact must not import the curator.
"""

from pathlib import Path, PurePosixPath
from types import SimpleNamespace

import lamindb as ln
import lamindb_setup
import pytest
from lamindb.core.storage.paths import AUTO_KEY_PREFIX
from lamindb.errors import FieldValidationError, InvalidArgument
from lamindb.models.artifact import get_relative_path_to_directory
from lamindb_setup.core.canonical_suffix import CanonicalSuffix

ln.settings.verbosity = "success"


@pytest.fixture
def tsv_file(tmp_path: Path) -> Path:
    filepath = tmp_path / "test.tsv"
    filepath.write_text("a\tb\n1\t2\n")
    return filepath


def test_basic_validation():
    # extra kwargs
    with pytest.raises(FieldValidationError):
        ln.Artifact("testpath.csv", description="test1b", extra_kwarg="extra")

    # > 1 args
    with pytest.raises(ValueError) as error:
        ln.Artifact("testpath.csv", "testpath.csv")
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed: path"

    # AUTO_KEY_PREFIX in key
    with pytest.raises(ValueError) as error:
        ln.Artifact(".gitignore", key=".lamindb/test_df.parquet")
    assert (
        error.exconly()
        == f"ValueError: Do not pass key that contains a managed storage path in `{AUTO_KEY_PREFIX}`"
    )

    # path that contains AUTO_KEY_PREFIX
    with pytest.raises(ValueError) as error:
        ln.Artifact(".lamindb/test_df.parquet", description="Test")
    assert (
        error.exconly()
        == f"ValueError: Do not pass path inside the `{AUTO_KEY_PREFIX}` directory for non-s3/gs paths."
    )


def test_cloud_path_init_missing_storage_raises(monkeypatch):
    path = "s3://missing-storage/.lamindb/test_df.parquet"
    monkeypatch.setattr(
        "lamindb.models.artifact.select_storage_or_parent", lambda _: None
    )
    with pytest.raises(ValueError) as error:
        ln.Artifact(path, description="test")
    assert (
        error.exconly()
        == f"ValueError: Registered storage location for path '{path}' not found."
    )


def test_cloud_path_init_missing_instance_slug_raises(monkeypatch):
    path = "s3://my-storage/.lamindb/test_df.parquet"
    storage_record = {"instance_uid": "missing-instance", "root": "s3://my-storage"}
    monkeypatch.setattr(
        "lamindb.models.artifact.select_storage_or_parent", lambda _: storage_record
    )
    monkeypatch.setattr(
        "lamindb.models.artifact.get_instance_slug_by_uid", lambda _: None
    )
    with pytest.raises(ValueError) as error:
        ln.Artifact(path, description="test")
    assert (
        error.exconly()
        == f"ValueError: Managing instance for storage location '{storage_record['root']}' not found."
    )


def test_cloud_path_init_missing_artifact_raises(monkeypatch):
    path = "s3://my-storage/.lamindb/test_df.parquet"
    monkeypatch.setattr(
        "lamindb.models.artifact.select_storage_or_parent",
        lambda _: {"instance_uid": "instance-uid", "root": "s3://my-storage"},
    )
    monkeypatch.setattr(
        "lamindb.models.artifact.get_instance_slug_by_uid",
        lambda _: "owner/instance",
    )

    class DummyConnection:
        def get(self, **_kwargs):
            raise ln.Artifact.DoesNotExist

    monkeypatch.setattr(
        ln.Artifact,
        "connect",
        classmethod(lambda cls, instance: DummyConnection()),
    )
    with pytest.raises(ValueError) as error:
        ln.Artifact(path, description="test")
    assert error.exconly() == f"ValueError: Artifact for path '{path}' not found."


_FOLDER_ARTIFACT_PATHS = [
    (
        "s3://my-storage/.lamindb/abcdefghijklmnop/subdir/file.txt",
        "s3://my-storage/.lamindb/abcdefghijklmnop",
    ),
    (
        "s3://my-storage/.lamindb/abcdefghijklmnop.zarr/group/arr",
        "s3://my-storage/.lamindb/abcdefghijklmnop.zarr",
    ),
]


@pytest.mark.parametrize("nested_path,folder_path", _FOLDER_ARTIFACT_PATHS)
def test_cloud_path_init_file_in_folder_artifact_current_instance(
    monkeypatch, nested_path, folder_path
):
    dummy = SimpleNamespace(uid="abcdefghijklmnop0000")
    looked_up = []

    def fake_get(cls, idlike=None, **kwargs):
        path = kwargs.get("path")
        looked_up.append(path)
        if path == folder_path:
            return dummy
        raise ln.Artifact.DoesNotExist

    inits = []

    monkeypatch.setattr(ln.Artifact, "get", classmethod(fake_get))
    monkeypatch.setattr(
        "lamindb.models.sqlrecord.init_self_from_db",
        lambda self, existing, db="default": inits.append(existing),
    )

    ln.Artifact(nested_path)
    assert looked_up == [nested_path, folder_path]
    assert inits == [dummy]


@pytest.mark.parametrize("nested_path,folder_path", _FOLDER_ARTIFACT_PATHS)
def test_cloud_path_init_file_in_folder_artifact_associated_instance(
    monkeypatch, nested_path, folder_path
):
    dummy = SimpleNamespace(
        uid="abcdefghijklmnop0000", _state=SimpleNamespace(db="owner/instance")
    )

    def fake_current_get(cls, idlike=None, **kwargs):
        raise ln.Artifact.DoesNotExist

    monkeypatch.setattr(ln.Artifact, "get", classmethod(fake_current_get))
    monkeypatch.setattr(
        "lamindb.models.artifact.select_storage_or_parent",
        lambda _: {"instance_uid": "instance-uid", "root": "s3://my-storage"},
    )
    monkeypatch.setattr(
        "lamindb.models.artifact.get_instance_slug_by_uid",
        lambda _: "owner/instance",
    )

    looked_up = []

    class DummyConnection:
        def get(self, **kwargs):
            path = kwargs.get("path")
            looked_up.append(path)
            if path == folder_path:
                return dummy
            raise ln.Artifact.DoesNotExist

    monkeypatch.setattr(
        ln.Artifact,
        "connect",
        classmethod(lambda cls, instance: DummyConnection()),
    )
    inits = []
    monkeypatch.setattr(
        "lamindb.models.sqlrecord.init_self_from_db",
        lambda self, existing, db="default": inits.append((existing, db)),
    )

    ln.Artifact(nested_path)
    assert looked_up == [nested_path, folder_path]
    assert inits == [(dummy, "owner/instance")]


def test_cloud_path_init_file_in_folder_artifact_missing_raises(monkeypatch):
    nested_path = "s3://my-storage/.lamindb/abcdefghijklmnop/subdir/file.txt"
    monkeypatch.setattr(
        "lamindb.models.artifact.select_storage_or_parent",
        lambda _: {"instance_uid": "instance-uid", "root": "s3://my-storage"},
    )
    monkeypatch.setattr(
        "lamindb.models.artifact.get_instance_slug_by_uid",
        lambda _: "owner/instance",
    )

    class DummyConnection:
        def get(self, **_kwargs):
            raise ln.Artifact.DoesNotExist

    monkeypatch.setattr(
        ln.Artifact,
        "connect",
        classmethod(lambda cls, instance: DummyConnection()),
    )
    with pytest.raises(ValueError) as error:
        ln.Artifact(nested_path, description="test")
    assert (
        error.exconly() == f"ValueError: Artifact for path '{nested_path}' not found."
    )


@pytest.mark.parametrize("key_is_virtual", [True, False])
@pytest.mark.parametrize("key", [None, "my_new_dir/my_artifact.csv", "nosuffix"])
@pytest.mark.parametrize("description", [None, "my description"])
def test_create_from_path_file(get_test_filepaths, key_is_virtual, key, description):
    ln.settings.creation._artifact_use_virtual_keys = key_is_virtual
    is_in_registered_storage = get_test_filepaths[0]
    root_dir = get_test_filepaths[1]
    test_filepath = get_test_filepaths[3]
    suffix = get_test_filepaths[4]  # path suffix
    if key is not None:
        key_suffix = CanonicalSuffix.from_path(PurePosixPath(key))
    else:
        key_suffix = None
    # this tests if insufficient information is being provided
    if key is None and not is_in_registered_storage and description is None:
        # this can fail because ln.track() might set a global run context
        # in that case, the Artifact would have a run that's not None and the
        # error below wouldn't be thrown
        with pytest.raises(ValueError) as error:
            artifact = ln.Artifact(test_filepath, key=key, description=description)
        assert (
            error.exconly()
            == "ValueError: Pass one of key, run or description as a parameter"
        )
        return None
    elif key is not None and suffix != key_suffix:
        with pytest.raises(InvalidArgument) as error:
            artifact = ln.Artifact(test_filepath, key=key, description=description)
        assert error.exconly() == (
            f"lamindb.errors.InvalidArgument: The passed path's suffix '{suffix}' must match the passed key's suffix '{key_suffix}'."
        )
        return None
    elif key is not None and is_in_registered_storage:
        inferred_key = get_relative_path_to_directory(
            path=test_filepath, directory=root_dir
        ).as_posix()
        try:
            artifact = ln.Artifact(test_filepath, key=key, description=description)
        except InvalidArgument as error:
            assert str(error) == (
                f"The path '{test_filepath}' is already in registered storage"
                f" '{root_dir.resolve().as_posix()}' with key '{inferred_key}'\nYou"
                f" passed conflicting key '{key}': please move the file before"
                " registering it."
            )
        return None
    else:
        artifact = ln.Artifact(test_filepath, key=key, description=description)
        assert artifact._state.adding  # make sure that this is a new file in the db
    assert (
        artifact.description is None
        if description is None
        else artifact.description == description
    )
    assert artifact.suffix == suffix
    assert artifact.n_files is None
    artifact.save()
    assert artifact.path.exists()
    # check get by path
    assert ln.Artifact.get(path=artifact.path) == artifact

    if key is None:
        assert (
            artifact.key == f"my_dir/my_file{suffix}"
            if is_in_registered_storage
            else artifact.key is None
        )
        if is_in_registered_storage:
            assert artifact.storage.root == root_dir.resolve().as_posix()
            assert artifact.path == test_filepath.resolve()
        else:
            assert artifact.storage.root == lamindb_setup.settings.storage.root_as_str
            assert (
                artifact.path
                == lamindb_setup.settings.storage.root
                / f".lamindb/{artifact.uid}{suffix}"
            )
    else:
        assert artifact.key == key
        assert artifact._key_is_virtual == key_is_virtual
        if is_in_registered_storage:
            # this would only hit if the key matches the correct key
            assert artifact.storage.root == root_dir.resolve().as_posix()
            assert (
                artifact.path == root_dir / f"{key}{suffix}" == test_filepath.resolve()
            )
        else:
            # file is moved into default storage
            if key_is_virtual:
                assert (
                    artifact.path
                    == lamindb_setup.settings.storage.root
                    / f".lamindb/{artifact.uid}{suffix}"
                )
            else:
                assert artifact.path == lamindb_setup.settings.storage.root / key
    # only delete from storage if a file copy took place
    delete_from_storage = str(test_filepath.resolve()) != str(artifact.path)
    artifact.delete(permanent=True, storage=delete_from_storage)
    ln.settings.creation._artifact_use_virtual_keys = True


@pytest.mark.parametrize("key_is_virtual", [True, False])
@pytest.mark.parametrize("key", [None, "my_new_file.tsv"])
def test_create_from_path_file_with_explicit_key_is_virtual(
    tsv_file, key_is_virtual, key
):
    artifact = ln.Artifact(
        tsv_file,
        description="test explicit key is virtual",
        key=key,
        key_is_virtual=key_is_virtual,
    )
    assert artifact.key == key
    assert artifact._key_is_virtual == key_is_virtual
    artifact.save()
    assert artifact.path.exists()

    root = lamindb_setup.settings.storage.root
    if not key_is_virtual and key is not None:
        assert artifact.path == root / key
    else:
        assert artifact.path == root / f".lamindb/{artifact.uid}.tsv"

    artifact.delete(permanent=True, storage=True)


def test_create_from_path_file_with_both_key_is_virtual_args_raises(tsv_file):
    with pytest.raises(ValueError) as error:
        ln.Artifact(
            tsv_file,
            key="my_new_file.tsv",
            key_is_virtual=True,
            _key_is_virtual=False,
        )
    assert (
        error.exconly()
        == "ValueError: Do not pass both key_is_virtual and _key_is_virtual."
    )


def test_create_from_empty_files_skips_hash_lookup(tmp_path):
    path_1 = tmp_path / "empty-1.txt"
    path_2 = tmp_path / "empty-2.txt"
    path_1.write_text("")
    path_2.write_text("")

    artifact_1 = ln.Artifact(path_1, key=f"{tmp_path.name}/empty-1.txt").save()
    artifact_2 = ln.Artifact(path_2, key=f"{tmp_path.name}/empty-2.txt")

    assert artifact_2.uid != artifact_1.uid
    assert artifact_2.key == f"{tmp_path.name}/empty-2.txt"
    assert artifact_2.hash == artifact_1.hash

    artifact_2.save()
    assert artifact_2.id != artifact_1.id

    artifact_2.delete(permanent=True)
    artifact_1.delete(permanent=True)


def test_existing_storage_path_skips_hash_lookup_by_default(tmp_path):
    storage_root = tmp_path / "registered-storage"
    storage_root.mkdir()
    filepath_1 = storage_root / "existing-1.txt"
    filepath_2 = storage_root / "existing-2.txt"
    filepath_1.write_text("same-content")
    filepath_2.write_text("same-content")
    storage = ln.Storage(root=storage_root.resolve().as_posix(), type="local").save()

    artifact_1 = ln.Artifact(filepath_1).save()
    artifact_2 = ln.Artifact(filepath_2, description="register sibling file")

    assert artifact_2._state.adding
    assert artifact_2.uid != artifact_1.uid
    assert artifact_2.hash == artifact_1.hash

    artifact_2.save()
    assert artifact_2.id != artifact_1.id
    assert artifact_2.key != artifact_1.key

    artifact_2.delete(permanent=True, storage=True)
    artifact_1.delete(permanent=True, storage=True)
    storage.delete()


def test_upload_checks_hash_by_default(tmp_path):
    filepath = tmp_path / "uploaded.txt"
    filepath.write_text("uploaded-content")

    artifact_1 = ln.Artifact(
        filepath, key="uploads/uploaded.txt", extra_data={"source": "manual"}
    ).save()
    artifact_2 = ln.Artifact(filepath)

    assert not artifact_2._state.adding
    assert artifact_2.id == artifact_1.id
    assert artifact_2.uid == artifact_1.uid
    assert artifact_2.extra_data == {"source": "manual"}

    artifact_1.delete(permanent=True)


def test_invalid_suffix_is_empty(tmp_path):
    filepath = tmp_path / "test.xyz"
    filepath.write_text("test-content")

    with pytest.raises(InvalidArgument) as error:
        ln.Artifact(filepath, key="uploads/test")
    assert error.exconly() == (
        "lamindb.errors.InvalidArgument: The passed path's suffix '.xyz' must match"
        " the passed key's suffix ''."
    )

    artifact = ln.Artifact(
        filepath, key="uploads/test.xyz", description="unknown text suffix"
    )
    artifact.save()

    assert artifact.suffix == ""
    assert artifact.key == "uploads/test.xyz"

    artifact.delete(permanent=True)


def test_existing_storage_can_force_hash_lookup(tmp_path):
    storage_root = tmp_path / "registered-storage-check"
    storage_root.mkdir()
    filepath = storage_root / "existing.txt"
    filepath.write_text("same-content")
    storage = ln.Storage(root=storage_root.resolve().as_posix(), type="local").save()

    artifact_1 = ln.Artifact(filepath).save()
    artifact_2 = ln.Artifact(filepath, skip_hash_lookup=False)

    assert not artifact_2._state.adding
    assert artifact_2.id == artifact_1.id

    artifact_1.delete(permanent=True, storage=False)
    storage.delete()


def test_skip_hash_lookup_true_on_upload_creates_new_artifact(tmp_path):
    filepath = tmp_path / "skip-true.txt"
    filepath.write_text("skip-true")

    # Use different keys: with the same key, save() resolves a (storage, key, hash)
    # collision back to the existing persisted artifact.
    artifact_1 = ln.Artifact(filepath, key="uploads/skip-true-a.txt").save()
    artifact_2 = ln.Artifact(
        filepath,
        key="uploads/skip-true-b.txt",
        skip_hash_lookup=True,
    )
    assert artifact_2._state.adding
    assert artifact_2.uid != artifact_1.uid
    artifact_2.save()
    assert artifact_2.id != artifact_1.id
    assert artifact_2.hash == artifact_1.hash
    assert artifact_2.key != artifact_1.key

    artifact_2.delete(permanent=True)
    artifact_1.delete(permanent=True)


def test_skip_hash_lookup_true_on_upload_same_key_resolves_collision_on_save(tmp_path):
    filepath = tmp_path / "skip-true-same-key.txt"
    filepath.write_text("skip-true")

    artifact_1 = ln.Artifact(filepath, key="uploads/skip-true-same-key.txt").save()
    artifact_2 = ln.Artifact(
        filepath,
        key="uploads/skip-true-same-key.txt",
        skip_hash_lookup=True,
    )
    # Constructor still creates a new unsaved object when hash lookup is skipped.
    assert artifact_2._state.adding
    assert artifact_2.uid != artifact_1.uid

    # On save(), (storage, key, hash) uniqueness resolves to the persisted record
    # in SQLRecord.save() (lamindb/models/sqlrecord.py), which calls
    # init_self_from_db(self, pre_existing_record) on hash/key collision.
    artifact_2.save()
    assert artifact_2.id == artifact_1.id
    assert artifact_2.uid == artifact_1.uid

    artifact_1.delete(permanent=True)


def test_stats_override_skips_hashing(tmp_path, monkeypatch):
    import lamindb.models.artifact as artifact_module

    filepath = tmp_path / "stats.txt"
    filepath.write_text("content that is not hashed")

    def no_hashing(*args, **kwargs):
        raise AssertionError("hash_file must not run when stats are passed")

    monkeypatch.setattr(artifact_module, "hash_file", no_hashing)
    stats = {"size": 3, "hash": "kAFQmDzST7DWlj99KOF_cg", "hash_type": "md5"}
    artifact = ln.Artifact(filepath, key="uploads/stats.txt", _stats=stats)
    assert artifact.size == 3
    assert artifact.hash == "kAFQmDzST7DWlj99KOF_cg"
    assert artifact._hash_type == "md5"
    assert artifact.n_files is None
    assert not artifact._overwrite_versions
    artifact.save()
    loaded = ln.Artifact.get(artifact.uid)
    assert loaded.hash == "kAFQmDzST7DWlj99KOF_cg"
    assert loaded.size == 3
    artifact.delete(permanent=True)


def test_stats_is_not_public(tmp_path):
    from lamindb.errors import FieldValidationError
    from lamindb.models.sqlrecord import _get_record_kwargs

    assert "stats" not in [name for name, _ in _get_record_kwargs(ln.Artifact)]
    assert "_stats" not in [name for name, _ in _get_record_kwargs(ln.Artifact)]
    filepath = tmp_path / "public.txt"
    filepath.write_text("public")
    with pytest.raises(FieldValidationError):
        ln.Artifact(filepath, key="uploads/public.txt", stats={"size": 6})


def test_stats_folder_sets_overwrite_versions(tmp_path):
    folder = tmp_path / "stats-folder"
    folder.mkdir()
    (folder / "a").write_text("a")
    stats = {
        "size": 2,
        "hash": "YNCwyFw7rBegqZLBVEckCg",
        "hash_type": "md5-d",
        "n_files": 2,
    }
    artifact = ln.Artifact(folder, key="uploads/stats-folder", _stats=stats)
    assert artifact.n_files == 2
    assert artifact._overwrite_versions


def test_stats_dedups_by_hash(tmp_path):
    filepath = tmp_path / "dedup.txt"
    filepath.write_text("dedup-content")
    artifact_1 = ln.Artifact(filepath, key="uploads/dedup.txt").save()

    other = tmp_path / "other.txt"
    other.write_text("different content, same declared hash")
    stats = {
        "size": artifact_1.size,
        "hash": artifact_1.hash,
        "hash_type": artifact_1._hash_type,
    }
    artifact_2 = ln.Artifact(
        other, key="uploads/other.txt", _stats=stats, skip_hash_lookup=False
    )
    assert not artifact_2._state.adding
    assert artifact_2.id == artifact_1.id

    artifact_1.delete(permanent=True)


def test_stats_without_hash_creates_new_artifact(tmp_path):
    filepath = tmp_path / "nohash.txt"
    filepath.write_text("nohash")
    artifact = ln.Artifact(filepath, key="uploads/nohash.txt", _stats={"size": 6})
    assert artifact._state.adding
    assert artifact.size == 6
    assert artifact.hash is None
    assert artifact._hash_type is None


def test_stats_invalid_raises(tmp_path):
    filepath = tmp_path / "bad.txt"
    filepath.write_text("bad")
    with pytest.raises(ValueError, match=r"_stats\['hash'\]"):
        ln.Artifact(
            filepath,
            key="uploads/bad.txt",
            _stats={"size": 3, "hash": "x" * 23, "hash_type": "md5"},
        )
    with pytest.raises(ValueError, match=r"_stats\['hash_type'\]"):
        ln.Artifact(
            filepath,
            key="uploads/bad.txt",
            _stats={"size": 3, "hash": "kAFQmDzST7DWlj99KOF_cg"},
        )
    with pytest.raises(ValueError, match=r"_stats\['size'\]"):
        ln.Artifact(filepath, key="uploads/bad.txt", _stats={"size": -1})
    with pytest.raises(ValueError, match="must be a dict"):
        ln.Artifact(filepath, key="uploads/bad.txt", _stats=["size"])
    with pytest.raises(ValueError, match="unknown keys"):
        ln.Artifact(filepath, key="uploads/bad.txt", _stats={"size": 3, "md5": "abc"})
    with pytest.raises(ValueError, match=r"requires _stats\['hash'\]"):
        ln.Artifact(filepath, key="uploads/bad.txt", _stats={"hash_type": "md5"})


@pytest.mark.parametrize("key", [None, "my_new_folder"])
def test_create_from_path_folder(get_test_filepaths, key):
    # get variables from fixture
    is_in_registered_storage = get_test_filepaths[0]
    test_dirpath = get_test_filepaths[2]
    hash_test_dir = get_test_filepaths[5]
    if key is None and not is_in_registered_storage:
        with pytest.raises(ValueError) as error:
            ln.Artifact(test_dirpath, key=key)
        assert error.exconly().startswith(
            "ValueError: Pass one of key, run or description as a parameter"
        )
        return None
    artifact1 = ln.Artifact(test_dirpath, key=key)
    if key is not None and is_in_registered_storage:
        assert artifact1._real_key is not None
        # should fail because we are passing a path in an existing storage with a virtual key
        with pytest.raises(ValueError) as error:
            ln.Artifact(test_dirpath, key=key, key_is_virtual=False)
        assert error.exconly().startswith(
            "ValueError: Passing a path in an existing storage with a virtual key and key_is_virtual=False is incompatible."
        )
    else:
        assert artifact1._real_key is None
    # check that passing key_is_virtual=True is incompatible with a path in an existing storage without a virtual key
    if key is None and is_in_registered_storage:
        with pytest.raises(ValueError) as error:
            ln.Artifact(test_dirpath, key=key, key_is_virtual=True)
        assert error.exconly().startswith(
            "ValueError: Passing a path in an existing storage without a virtual key and key_is_virtual=True is incompatible."
        )
    assert artifact1.n_files == 3
    assert artifact1.hash == hash_test_dir
    assert artifact1._state.adding
    assert artifact1.description is None
    assert artifact1.path.exists()
    artifact1.save()

    # run tests on re-creating the Artifact
    artifact2 = ln.Artifact(test_dirpath, key=key, description="something")
    if is_in_registered_storage:
        assert artifact2._state.adding
        assert artifact1.uid != artifact2.uid
    else:
        assert not artifact2._state.adding
        assert artifact1.id == artifact2.id
        assert artifact1.uid == artifact2.uid
    assert artifact1.storage == artifact2.storage
    assert artifact2.path.exists()
    assert artifact2.description == "something"
    artifact2.save()
    if is_in_registered_storage:
        # init skips hash lookup by default, but save() still resolves same
        # storage+key+hash to the persisted record.
        assert artifact1.id == artifact2.id

    # now put another file in the test directory

    # create a first file
    test_filepath_added = test_dirpath / "my_file_added.txt"
    test_filepath_added.write_text("2")
    artifact3 = ln.Artifact(test_dirpath, key=key, revises=artifact1)
    assert artifact3.n_files == 4
    assert artifact3.hash != hash_test_dir
    assert artifact3._state.adding
    assert artifact3.description is None
    assert artifact3.path.exists()
    artifact3.save()

    # the state of artifact1 is lost, because artifact3 is stored at the same path
    assert artifact3._single_path_across_versions
    assert artifact1._single_path_across_versions
    assert artifact3.path == artifact1.path
    test_filepath_added.unlink()

    # delete the artifact
    artifact2.delete(permanent=True, storage=False)
    artifact3.delete(permanent=True, storage=False)


def test_create_from_path_overwrite_versions_false(get_test_filepaths):
    # get variables from fixture
    is_in_registered_storage = get_test_filepaths[0]
    test_dirpath = get_test_filepaths[2]
    hash_test_dir = get_test_filepaths[5]
    if is_in_registered_storage:
        return
    with pytest.warns(DeprecationWarning, match="overwrite_versions"):
        artifact1 = ln.Artifact(
            test_dirpath, key="my_folder", overwrite_versions=False
        ).save()
    assert artifact1.hash == hash_test_dir
    assert artifact1._single_path_across_versions is False
    with pytest.warns(DeprecationWarning, match="_single_path_across_versions"):
        assert artifact1.overwrite_versions is False
    # skip artifact2 because we already test this above
    # create a first file
    test_filepath_added = test_dirpath / "my_file_added.txt"
    test_filepath_added.write_text("2")
    with pytest.warns(DeprecationWarning, match="overwrite_versions"):
        artifact3 = ln.Artifact(test_dirpath, key="my_folder", overwrite_versions=False)
    assert artifact3.hash != hash_test_dir
    artifact3.save()
    # the state of artifact1 is lost, because artifact3 is stored at the same path
    assert not artifact3._single_path_across_versions
    assert not artifact1._single_path_across_versions
    assert artifact3.path != artifact1.path
    test_filepath_added.unlink()
    artifact1.delete(permanent=True, storage=False)
    artifact3.delete(permanent=True, storage=False)


def test_create_from_path_set_branch():
    branch = ln.Branch(name="contrib1").save()
    artifact1 = ln.Artifact(".gitignore", key="test", branch=branch).save()
    # check hash lookup on different branch
    artifact2 = ln.Artifact(".gitignore", key="test1")
    assert artifact1 == artifact2
    # cleanup
    artifact1.delete(permanent=True)
    branch.delete(permanent=True)


@pytest.mark.parametrize("key", [None, "my_new_folder"])
def test_from_dir(get_test_filepaths, key):
    is_in_registered_storage = get_test_filepaths[0]
    test_dirpath = get_test_filepaths[2]
    # the directory contains 3 files, two of them are duplicated
    artifacts = ln.Artifact.from_dir(test_dirpath, key=key)
    for artifact in artifacts:
        if key is not None and is_in_registered_storage:
            assert artifact._real_key is not None
        else:
            assert artifact._real_key is None
    # we only return the duplicated ones
    hashes = [artifact.hash for artifact in artifacts if artifact.hash is not None]
    uids = [artifact.uid for artifact in artifacts]
    assert len(set(hashes)) == len(hashes)
    ln.UPath(test_dirpath).view_tree()
    # now save
    artifacts.save()
    # now run again; in existing storage this should skip hash lookup by default
    artifacts = ln.Artifact.from_dir(test_dirpath, key=key)
    assert len(artifacts) == 2
    assert len({artifact.uid for artifact in artifacts}) == len(hashes)
    if is_in_registered_storage:
        assert all(artifact._state.adding for artifact in artifacts)
    else:
        assert all(not artifact._state.adding for artifact in artifacts)
    queried_artifacts = ln.Artifact.filter(uid__in=uids)
    for artifact in queried_artifacts:
        artifact.delete(permanent=True, storage=False)
