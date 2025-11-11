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


def test_dataframe_suffix(example_dataframe: pd.DataFrame):
    df = example_dataframe
    artifact = ln.Artifact.from_dataframe(df, key="test_.parquet")
    assert artifact.suffix == ".parquet"

    with pytest.raises(ln.errors.InvalidArgument) as error:
        artifact = ln.Artifact.from_dataframe(df, key="test_.def")
    assert (
        error.exconly().partition(",")[0]
        == "lamindb.errors.InvalidArgument: The passed key's suffix '.def' must match the passed path's suffix '.parquet'."
    )


def test_artifact_dataframe_with_features(example_dataframe: pd.DataFrame):
    """Test column names encoding when features with the same names are present."""
    artifact = ln.Artifact.from_dataframe(example_dataframe, key="df.parquet").save()
    id_feature = ln.Feature(name="id", dtype="int").save()
    uid_feature = ln.Feature(name="uid", dtype="str").save()
    artifact.features.add_values({"id": 1, "uid": "test-uid"})
    df = ln.Artifact.filter(key="df.parquet").to_dataframe(
        include=["description"], features=True
    )
    assert df.index.name == "__lamindb_artifact_id__"
    assert df.columns.tolist() == [
        "__lamindb_artifact_uid__",
        "key",
        "id",
        "uid",
        "description",
    ]
    assert df.iloc[0]["id"] == 1
    assert df.iloc[0]["uid"] == "test-uid"

    artifact.delete(permanent=True)
    id_feature.delete(permanent=True)
    uid_feature.delete(permanent=True)


def test_artifact_from_dataframe_with_schema(example_dataframe: pd.DataFrame):
    df = example_dataframe
    feat1 = ln.Feature(name="feat1", dtype=int).save()
    artifact = ln.Artifact.from_dataframe(
        df, key="test_df.parquet", schema="valid_features"
    ).save()
    assert artifact.schema == ln.examples.schemas.valid_features()
    inferred_schema_link = artifact.feature_sets.through.get(artifact_id=artifact.id)
    assert inferred_schema_link.slot == "columns"
    assert inferred_schema_link.schema.members.count() == 1
    assert inferred_schema_link.schema.members.first() == feat1
    inferred_schema = inferred_schema_link.schema
    inferred_schema_link.delete()
    inferred_schema.delete(permanent=True)
    feat1.delete(permanent=True)
    artifact.delete(permanent=True)
