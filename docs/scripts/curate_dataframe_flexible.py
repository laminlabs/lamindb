import lamindb as ln

ln.core.datasets.mini_immuno.define_features_labels()
df = ln.core.datasets.small_dataset1(otype="DataFrame")
schema = ln.schemas.simple.valid_features()
artifact = ln.Artifact.from_df(
    df, key="examples/dataset1.parquet", schema=schema
).save()
artifact.describe()
