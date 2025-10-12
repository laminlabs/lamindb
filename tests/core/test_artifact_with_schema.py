# Note: Almost all logic for schema-based validation is handled in the curators test suite
# This here only covers external feature annotation and validation
from pathlib import Path

import lamindb as ln
import pandas as pd
import pytest


@pytest.fixture(scope="module")
def two_features():
    feat1 = ln.Feature(name="species", dtype=str).save()
    feat2 = ln.Feature(name="split", dtype=str).save()
    yield feat1, feat2
    feat1.delete(permanent=True)
    feat2.delete(permanent=True)


@pytest.mark.parametrize("use_schema", [True, False])
def test_create_external_schema(
    tsv_file: Path, use_schema: bool, two_features: tuple[ln.Feature, ln.Feature]
):
    feat1, feat2 = two_features
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
    artifact.delete(permanent=True)
    if use_schema:
        schema.delete(permanent=True)


def test_from_dataframe_external_schema(
    df: pd.DataFrame, two_features: tuple[ln.Feature, ln.Feature]
):
    feat1, feat2 = two_features
    external_schema = ln.Schema(itype=ln.Feature).save()
    schema = ln.Schema(
        features=[feat1, feat2],
        slots={"__external__": external_schema},
        otype="DataFrame",
    ).save()
    artifact = ln.Artifact.from_dataframe(
        df,
        key="test_df_with_external_features.parquet",
        features={"species": "bird", "split": "train"},
        schema=schema,
    ).save()
    assert artifact.features.get_values() == {"species": "bird", "split": "train"}
    artifact.delete(permanent=True)
    schema.delete(permanent=True)
    external_schema.delete(permanent=True)
