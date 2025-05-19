import lamindb as ln

ln.core.datasets.mini_immuno.define_features_labels()
schema = ln.examples.schemas.valid_features()
df = ln.core.datasets.small_dataset1(otype="DataFrame")
artifact = ln.Artifact.from_df(
    df, key="examples/dataset1.parquet", schema=schema
).save()
artifact.describe()
