import lamindb as ln

from .define_schema_df_metadata import study_metadata_schema

df = ln.examples.datasets.mini_immuno.get_dataset1(otype="DataFrame")
schema = ln.Schema(
    features=[ln.Feature(name="perturbation", dtype="str").save()],
    slots={"attrs": study_metadata_schema},
    otype="DataFrame",
).save()
curator = ln.curators.DataFrameCurator(df, schema=schema)
curator.validate()
artifact = curator.save_artifact(key="examples/df_with_attrs.parquet")
artifact.describe()
