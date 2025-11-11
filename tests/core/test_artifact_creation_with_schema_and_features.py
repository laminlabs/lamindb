# Note: Almost all logic for schema-based validation is handled in the curators test suite
# This here only covers external feature annotation and validation

import lamindb as ln
import pandas as pd
import pytest


@pytest.fixture(scope="module")
def two_internal_features():
    feat1 = ln.Feature(name="feat1", dtype=int).save()
    feat2 = ln.Feature(name="feat2", dtype=int).save()
    yield feat1, feat2
    feat1.delete(permanent=True)
    feat2.delete(permanent=True)


@pytest.fixture(scope="module")
def two_external_features():
    feature_a = ln.Feature(name="feature_a", dtype=str).save()
    feature_b = ln.Feature(name="feature_b", dtype=str).save()
    yield feature_a, feature_b
    feature_a.delete(permanent=True)
    feature_b.delete(permanent=True)


@pytest.mark.parametrize("use_schema", [True, False])
def test_create_artifact_with_external_feature_annotations(
    use_schema: bool,
    two_external_features: tuple[ln.Feature, ln.Feature],
):
    feat1, feat2 = two_external_features
    if use_schema:
        schema = ln.Schema(features=[feat1, feat2]).save()
    else:
        schema = None
    artifact = ln.Artifact(
        ".gitignore",
        key="test_file",
        features={"feature_a": "x", "feature_b": "y"},
        schema=schema,
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}
    assert artifact.schema == schema
    # repeat to check idempotency (requires set_values() instead of add_values())
    artifact = ln.Artifact(
        ".gitignore",
        key="test_file",
        features={"feature_a": "x", "feature_b": "y"},
        schema=schema,
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}
    assert artifact.schema == schema
    if use_schema:
        with pytest.raises(ValueError) as error:
            artifact.features.remove_values("feature_a", value="x")
        assert (
            "Cannot remove values if artifact has external schema." in error.exconly()
        )
    else:
        artifact.features.remove_values("feature_a", value="x")
        assert artifact.features.get_values() == {"feature_b": "y"}
    artifact.delete(permanent=True)
    if use_schema:
        schema.delete(permanent=True)


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


def test_from_dataframe_with_external_schema(
    example_dataframe: pd.DataFrame,
    two_external_features: tuple[ln.Feature, ln.Feature],
    two_internal_features: tuple[ln.Feature, ln.Feature],
):
    df = example_dataframe
    feat1, feat2 = two_internal_features
    featA, featB = two_external_features
    schema_external = ln.Schema(features=[featA, featB]).save()

    # Case 1: wrong internal features for this dataframe
    schema_with_mistake = ln.Schema(
        features=[featA, featB],
        slots={"__external__": schema_external},
        otype="DataFrame",
    ).save()
    with pytest.raises(ln.errors.ValidationError) as error:
        artifact = ln.Artifact.from_dataframe(
            df,
            key="test_df_with_external_features.parquet",
            features={"feature_a": "x", "feature_b": "y"},
            schema=schema_with_mistake,
        ).save()
    assert "COLUMN_NOT_IN_DATAFRAME" in error.exconly()

    # alternative via DataFrameCurator directly
    with pytest.raises(ln.errors.ValidationError) as error:
        ln.curators.DataFrameCurator(
            df,
            schema=schema_with_mistake,
        ).validate()
    assert "COLUMN_NOT_IN_DATAFRAME" in error.exconly()

    # Case 2: no schema for external features provided
    schema_no_external = ln.Schema(features=[feat1, feat2]).save()
    artifact = ln.Artifact.from_dataframe(
        df,
        key="test_df_with_external_features.parquet",
        features={"feature_a": "x", "feature_b": "y"},
        schema=schema_no_external,
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}
    artifact.delete(permanent=True)

    # alternative via DataFrameCurator directly
    curator = ln.curators.DataFrameCurator(
        df,
        schema=schema_no_external,
        features={"feature_a": "x", "feature_b": "y"},
    )
    artifact = curator.save_artifact(
        key="test_df_with_external_features.parquet",
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}
    artifact.delete(permanent=True)

    # Case 3: correct external schema
    schema_correct_external = ln.Schema(
        features=[feat1, feat2],
        slots={"__external__": schema_external},
        otype="DataFrame",
    ).save()
    artifact = ln.Artifact.from_dataframe(
        df,
        key="test_df_with_external_features.parquet",
        features={"feature_a": "x", "feature_b": "y"},
        schema=schema_correct_external,
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}
    assert (
        artifact.features.describe(return_str=True)
        == """\
Artifact: test_df_with_external_features.parquet (0000)
├── Dataset features
│   └── columns (2)
│       feat1               int
│       feat2               int
└── External features
    └── feature_a           str                     x
        feature_b           str                     y"""
    )
    with pytest.raises(ValueError) as error:
        artifact.features.remove_values("feature_a", value="x")
    assert "Cannot remove values if artifact has external schema." in error.exconly()
    artifact.delete(permanent=True)

    # alternative via DataFrameCurator directly
    curator = ln.curators.DataFrameCurator(
        df,
        schema=schema_correct_external,
        features={"feature_a": "x", "feature_b": "y"},
    )
    artifact = curator.save_artifact(
        key="test_df_with_external_features.parquet",
    ).save()
    assert artifact.features.get_values() == {"feature_a": "x", "feature_b": "y"}

    # clean up everything
    inferred_schema = artifact.feature_sets.all()[0]
    artifact.feature_sets.remove(inferred_schema.id)
    inferred_schema.delete(permanent=True)
    artifact.delete(permanent=True)
    schema_with_mistake.delete(permanent=True)
    schema_no_external.delete(permanent=True)
    schema_correct_external.delete(permanent=True)
    schema_external.delete(permanent=True)
