import lamindb as ln

schema = ln.examples.datasets.mini_immuno.define_mini_immuno_schema_flexible()
df = ln.examples.datasets.mini_immuno.get_dataset1(otype="DataFrame")
df.pop("donor")  # remove donor column to trigger validation error
try:
    artifact = ln.Artifact.from_dataframe(
        df, key="examples/dataset1.parquet", schema=schema
    ).save()
except ln.errors.ValidationError as error:
    print(error)
