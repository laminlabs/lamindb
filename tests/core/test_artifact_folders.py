import lamindb as ln
import pytest
from lamindb.errors import InvalidArgument


@pytest.mark.parametrize("key", [None, "my_new_folder"])
def test_folder_like_artifact(get_test_filepaths, key):
    # get variables from fixture
    is_in_registered_storage = get_test_filepaths[0]
    test_dirpath = get_test_filepaths[2]
    hash_test_dir = get_test_filepaths[5]

    # run tests on initial Artifact creation
    if key is not None and is_in_registered_storage:
        with pytest.raises(InvalidArgument) as error:
            ln.Artifact(test_dirpath, key=key)
        assert error.exconly().startswith(
            "lamindb.errors.InvalidArgument: The path"  # The path {data} is already in registered storage
        )
        return None
    if key is None and not is_in_registered_storage:
        with pytest.raises(ValueError) as error:
            ln.Artifact(test_dirpath, key=key)
        assert error.exconly().startswith(
            "ValueError: Pass one of key, run or description as a parameter"
        )
        return None
    artifact1 = ln.Artifact(test_dirpath, key=key)
    assert artifact1.n_files == 3
    assert artifact1.hash == hash_test_dir
    assert artifact1._state.adding
    assert artifact1.description is None
    assert artifact1.path.exists()
    artifact1.save()

    # run tests on re-creating the Artifact
    artifact2 = ln.Artifact(test_dirpath, key=key, description="something")
    assert not artifact2._state.adding
    assert artifact1.id == artifact2.id
    assert artifact1.uid == artifact2.uid
    assert artifact1.storage == artifact2.storage
    assert artifact2.path.exists()
    assert artifact2.description == "something"

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
    assert artifact3.path == artifact1.path
    test_filepath_added.unlink()

    # delete the artifact
    artifact2.delete(permanent=True, storage=False)
    artifact3.delete(permanent=True, storage=False)
