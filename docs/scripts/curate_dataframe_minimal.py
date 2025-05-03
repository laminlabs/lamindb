import lamindb as ln

schema = ln.core.datasets.mini_immuno.define_mini_immuno_schema_flexible()
df = ln.core.datasets.small_dataset1(otype="DataFrame")
curator = ln.curators.DataFrameCurator(df, schema)
artifact = curator.save_artifact(key="examples/dataset1.parquet")
assert artifact.schema == schema
