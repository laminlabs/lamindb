import lamindb as ln

schema = ln.schemas.dataframe.valid_features()

# curate a DataFrame
df = ln.core.datasets.small_dataset1(otype="DataFrame")
curator = ln.curators.DataFrameCurator(df, schema)
artifact = curator.save_artifact(key="example_datasets/dataset1.parquet")
assert artifact.schema == schema
