import lamindb as ln
import pandas as pd
import pytest


def test_create_from_dataframe(example_dataframe: pd.DataFrame):
    df = example_dataframe
    artifact = ln.Artifact.from_dataframe(df, description="test1")
    assert artifact.description == "test1"
    assert artifact.key is None
    assert artifact.otype == "DataFrame"
    assert artifact.kind == "dataset"
    assert artifact.n_observations == 2
    assert hasattr(artifact, "_local_filepath")
    artifact.key = "my-test-dataset"  # try changing key
    with pytest.raises(ln.errors.InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: The suffix '' of the provided key is incorrect, it should be '.parquet'."
    )
    artifact.key = None  # restore
    artifact.suffix = ".whatever"  # try changing suffix
    with pytest.raises(ln.errors.InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: Changing the `.suffix` of an artifact is not allowed! You tried to change it from '.parquet' to '.whatever'."
    )
    artifact.suffix = ".parquet"
    artifact.save()
    # check that the local filepath has been cleared
    assert not hasattr(artifact, "_local_filepath")
    del artifact

    # now get an artifact from the database
    artifact = ln.Artifact.get(description="test1")

    artifact.suffix = ".whatever"  # try changing suffix
    with pytest.raises(ln.errors.InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: Changing the `.suffix` of an artifact is not allowed! You tried to change it from '.parquet' to '.whatever'."
    )
    artifact.suffix = ".parquet"

    # coming from `key is None` that setting a key with different suffix is not allowed
    artifact.key = "my-test-dataset.suffix"
    with pytest.raises(ln.errors.InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: The suffix '.suffix' of the provided key is incorrect, it should be '.parquet'."
    )

    # coming from `key is None` test with no suffix
    artifact.key = "my-test-dataset"
    with pytest.raises(ln.errors.InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: The suffix '' of the provided key is incorrect, it should be '.parquet'."
    )

    # try a joint update where both suffix and key are changed, this previously enabled to create a corrupted state
    artifact.key = "my-test-dataset"
    artifact.suffix = ""
    with pytest.raises(ln.errors.InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: Changing the `.suffix` of an artifact is not allowed! You tried to change it from '.parquet' to ''."
    )

    # because this is a parquet artifact, we can set a key with a .parquet suffix
    artifact.suffix = ".parquet"  # restore proper suffix
    artifact.key = "my-test-dataset.parquet"
    artifact.save()
    assert artifact.key == "my-test-dataset.parquet"

    # coming from a .parquet key, test changing the key to no suffix
    artifact.key = "my-test-dataset"
    with pytest.raises(ln.errors.InvalidArgument) as error:
        artifact.save()
    assert (
        error.exconly()
        == "lamindb.errors.InvalidArgument: The suffix '' of the provided key is incorrect, it should be '.parquet'."
    )

    artifact.delete(permanent=True)


def test_revise_recreate_artifact(example_dataframe: pd.DataFrame, ccaplog):
    df = example_dataframe
    # attempt to create a file with an invalid version
    with pytest.raises(ValueError) as error:
        artifact = ln.Artifact.from_dataframe(df, description="test", version=0)
    assert (
        error.exconly()
        == "ValueError: `version` parameter must be `None` or `str`, e.g., '0.1', '1',"
        " '2', etc."
    )

    # create a file and tag it with a version
    key = "my-test-dataset.parquet"
    artifact = ln.Artifact.from_dataframe(df, key=key, description="test", version="1")
    assert artifact.version == "1"
    assert artifact.uid.endswith("0000")
    assert artifact.path.exists()  # because of cache file already exists
    artifact.save()
    assert artifact.path.exists()
    assert artifact.suffix == ".parquet"

    with pytest.raises(ValueError) as error:
        artifact_r2 = ln.Artifact.from_dataframe(df, revises=artifact, version="1")
    assert (
        error.exconly()
        == "ValueError: Please change the version tag or leave it `None`, '1' is already taken"
    )

    # create new file from old file
    df.iloc[0, 0] = 99  # mutate dataframe so that hash lookup doesn't trigger
    artifact_r2 = ln.Artifact.from_dataframe(df, revises=artifact)
    assert artifact_r2.stem_uid == artifact.stem_uid
    assert artifact_r2.uid.endswith("0001")
    # call this again
    artifact_r2 = ln.Artifact.from_dataframe(df, revises=artifact)
    assert artifact_r2.uid.endswith("0001")
    assert artifact_r2.stem_uid == artifact.stem_uid
    assert artifact_r2.version is None
    assert artifact_r2.key == key
    assert artifact.suffix == ".parquet"
    assert artifact_r2.description == "test"
    assert artifact_r2._revises is not None
    artifact_r2.save()
    assert artifact_r2.path.exists()
    assert artifact_r2._revises is None

    # create new file from newly versioned file
    df.iloc[0, 0] = 0  # mutate dataframe so that hash lookup doesn't trigger
    artifact_r3 = ln.Artifact.from_dataframe(
        df, description="test1", revises=artifact_r2, version="2"
    )
    assert artifact_r3.uid.endswith("0002")
    assert artifact_r3.stem_uid == artifact.stem_uid
    assert artifact_r3.version == "2"
    assert artifact_r3.description == "test1"

    # revise by matching on `key`
    df.iloc[0, 0] = 100  # mutate dataframe so that hash lookup doesn't trigger
    artifact_r3 = ln.Artifact.from_dataframe(
        df, description="test1", key=key, version="2"
    )
    assert artifact_r3.uid.endswith("0002")
    assert artifact_r3.stem_uid == artifact.stem_uid
    assert artifact_r3.key == key
    assert artifact_r3.version == "2"
    assert artifact_r3.description == "test1"
    assert artifact_r3.is_latest
    assert artifact_r2.is_latest
    artifact_r3.save()
    # now r2 is no longer the latest version, but need to re-fresh from db
    artifact_r2 = ln.Artifact.get(artifact_r2.uid)
    assert not artifact_r2.is_latest

    # re-create based on hash when artifact_r3 is in trash
    artifact_r3.delete()
    artifact_new = ln.Artifact.from_dataframe(
        df,
        key="my-test-dataset1.parquet",
    )
    assert artifact_new != artifact_r3
    assert artifact_new.hash == artifact_r3.hash
    assert artifact_new.key == "my-test-dataset1.parquet"
    artifact_r3.restore()  # restore from trash

    # re-create based on hash while providing a different key
    artifact_new = ln.Artifact.from_dataframe(
        df,
        key="my-test-dataset1.parquet",
        description="test1 updated",
    )
    assert artifact_new == artifact_r3
    assert artifact_new.hash == artifact_r3.hash
    assert artifact_new.key == key  # old key
    assert artifact_new.description == "test1 updated"

    # re-create while skipping hash lookup with different key
    artifact_v4 = ln.Artifact.from_dataframe(
        df,
        key="my-test-dataset1.parquet",
        skip_hash_lookup=True,
    )
    assert artifact_v4.uid != artifact_r3.uid
    assert artifact_v4.hash == artifact_r3.hash
    assert artifact_v4.key == "my-test-dataset1.parquet"
    artifact_v4.save()  # this just saves a duplicated file

    # re-create while skipping hash lookup with same key
    artifact_new = ln.Artifact.from_dataframe(
        df,
        key="my-test-dataset1.parquet",
        skip_hash_lookup=True,
    )
    assert artifact_new.uid != artifact_v4.uid
    assert artifact_new.stem_uid == artifact_v4.stem_uid
    assert artifact_new.hash == artifact_v4.hash
    artifact_new.save()  # should now violate unique constraint, falls back artifact_v4
    assert artifact_new.uid == artifact_v4.uid

    # re-create while skipping hash lookup artifact, move to trash before
    artifact_v4.delete()
    artifact_new = ln.Artifact.from_dataframe(
        df,
        key="my-test-dataset1.parquet",
        skip_hash_lookup=True,
    )
    assert artifact_new.uid != artifact_v4.uid
    assert artifact_new.key == "my-test-dataset1.parquet"
    assert "returning artifact from trash" not in ccaplog.text
    artifact_new.save()  # should now violate unique constraint, retrieve artifact_v4 from trash
    assert "returning artifact from trash" in ccaplog.text
    assert artifact_new.uid == artifact_v4.uid
    assert artifact_new.branch_id == 1  # restored to default branch

    with pytest.raises(TypeError) as error:
        ln.Artifact.from_dataframe(
            df, description="test1a", revises=ln.Record(name="test")
        )
    assert error.exconly() == "TypeError: `revises` has to be of type `Artifact`"

    artifact_r3.delete(permanent=True)
    artifact_r2.delete(permanent=True)
    artifact.delete(permanent=True)

    # unversioned file
    artifact = ln.Artifact.from_dataframe(df, description="test2")
    assert artifact.version is None

    # what happens if we don't save the old file?
    # add a test for it!
    artifact.save()

    # create new file from old file
    df.iloc[0, 0] = 101  # mutate dataframe so that hash lookup doesn't trigger
    new_artifact = ln.Artifact.from_dataframe(df, revises=artifact)
    assert artifact.version is None
    assert new_artifact.stem_uid == artifact.stem_uid
    assert new_artifact.version is None
    assert new_artifact.description == artifact.description

    artifact.delete()

    artifact_from_trash = ln.Artifact.get(artifact.uid[:-4])  # query with stem uid
    assert artifact_from_trash.branch_id == -1

    artifact.delete(permanent=True)  # permanent deletion
