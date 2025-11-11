import lamindb as ln
import pandas as pd


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
