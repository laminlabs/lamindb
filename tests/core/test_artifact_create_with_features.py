# Note: Almost all logic for schema-based validation is handled in the curators test suite
# This here only covers external feature annotation and validation
from pathlib import Path

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
    featA = ln.Feature(name="species", dtype=str).save()
    featB = ln.Feature(name="split", dtype=str).save()
    yield featA, featB
    featA.delete(permanent=True)
    featB.delete(permanent=True)


@pytest.mark.parametrize("use_schema", [True, False])
def test_create_artifact_with_external_feature_annotations(
    tsv_file: Path,
    use_schema: bool,
    two_external_features: tuple[ln.Feature, ln.Feature],
):
    feat1, feat2 = two_external_features
    if use_schema:
        schema = ln.Schema(features=[feat1, feat2]).save()
    else:
        schema = None
    artifact = ln.Artifact(
        tsv_file,
        key="test.tsv",
        features={"species": "bird", "split": "train"},
        schema=schema,
    ).save()
    assert artifact.features.get_values() == {"species": "bird", "split": "train"}
    assert artifact.schema == schema
    artifact.delete(permanent=True)
    if use_schema:
        schema.delete(permanent=True)


def test_from_dataframe_with_external_features_and_schema(
    df: pd.DataFrame,
    two_external_features: tuple[ln.Feature, ln.Feature],
    two_internal_features: tuple[ln.Feature, ln.Feature],
):
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
            features={"species": "bird", "split": "train"},
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
        features={"species": "bird", "split": "train"},
        schema=schema_no_external,
    ).save()
    assert artifact.features.get_values() == {"species": "bird", "split": "train"}
    artifact.delete(permanent=True)

    # alternative via DataFrameCurator directly
    curator = ln.curators.DataFrameCurator(
        df,
        schema=schema_no_external,
        features={"species": "bird", "split": "train"},
    )
    artifact = curator.save_artifact(
        key="test_df_with_external_features.parquet",
    ).save()
    assert artifact.features.get_values() == {"species": "bird", "split": "train"}
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
        features={"species": "bird", "split": "train"},
        schema=schema_correct_external,
    ).save()
    assert artifact.features.get_values() == {"species": "bird", "split": "train"}
    artifact.delete(permanent=True)

    # alternative via DataFrameCurator directly
    curator = ln.curators.DataFrameCurator(
        df,
        schema=schema_correct_external,
        features={"species": "bird", "split": "train"},
    )
    artifact = curator.save_artifact(
        key="test_df_with_external_features.parquet",
    ).save()
    assert artifact.features.get_values() == {"species": "bird", "split": "train"}

    # clean up everything
    inferred_schema = artifact.feature_sets.all()[0]
    artifact.feature_sets.remove(inferred_schema.id)
    inferred_schema.delete(permanent=True)
    artifact.delete(permanent=True)
    schema_with_mistake.delete(permanent=True)
    schema_no_external.delete(permanent=True)
    schema_correct_external.delete(permanent=True)
    schema_external.delete(permanent=True)
