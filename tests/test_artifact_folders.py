import lamindb as ln
import pytest


@pytest.mark.parametrize("key", [None, "my_new_folder"])
def test_folder_like_artifact(get_test_filepaths, key):
    isin_existing_storage = get_test_filepaths[0]
    test_dirpath = get_test_filepaths[2]
    hash_test_dir = get_test_filepaths[5]
    if key is not None and isin_existing_storage:
        with pytest.raises(ValueError) as error:
            ln.Artifact(test_dirpath, key=key)
        assert error.exconly().startswith(
            "ValueError: The path"  # The path {data} is already in registered storage
        )
        return None
    if key is None and not isin_existing_storage:
        with pytest.raises(ValueError) as error:
            ln.Artifact(test_dirpath, key=key)
        assert error.exconly().startswith(
            "ValueError: Pass one of key, run or description as a parameter"
        )
        return None
    artifact = ln.Artifact(test_dirpath, key=key)
    assert artifact.n_objects == 3
    assert artifact.hash == hash_test_dir
    assert artifact._state.adding
    assert artifact.description is None
    assert artifact.path.exists()
    artifact.save()
    # now run again, because now we'll have hash-based lookup!
    artifact2 = ln.Artifact(test_dirpath, key=key, description="something")
    assert not artifact2._state.adding
    assert artifact.id == artifact2.id
    assert artifact.uid == artifact2.uid
    assert artifact.storage == artifact2.storage
    assert artifact2.path.exists()
    assert artifact2.description == "something"
    artifact2.delete(permanent=True, storage=False)


# also see test in lamindb-setup/tests/storage/test_storage_stats.py
# there is also a test for GCP there
def test_folder_like_artifact_s3():
    study0_data = ln.Artifact(
        "s3://lamindb-dev-datasets/iris_studies/study0_raw_images"
    )
    assert study0_data.hash == "wVYKPpEsmmrqSpAZIRXCFg"
    assert study0_data.hash_type == "md5-d"
    assert study0_data.n_objects == 51
