import lamindb as ln

ln.examples.datasets.mini_immuno.define_features_labels()
df = ln.examples.datasets.mini_immuno.get_dataset1(otype="DataFrame")
artifact = ln.Artifact.from_dataframe(
    df, key="examples/dataset1.parquet", schema="valid_features"
).save()
artifact.describe()
