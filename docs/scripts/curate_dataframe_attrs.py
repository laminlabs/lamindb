import lamindb as ln

from .define_unstructured_schema import study_metadata_schema

df = ln.core.datasets.mini_immuno.get_dataset1(otype="DataFrame")
df.attrs = {"temperature": 21.6, "experiment_id": "EXP001"}

schema = ln.Schema(
    features=[ln.Feature(name="perturbation", dtype="str").save()],
    slots={"attrs": study_metadata_schema},
    otype="DataFrame",
).save()

curator = ln.curators.DataFrameCurator(df, schema=schema)
curator.validate()
artifact = curator.save_artifact(key="examples/df_with_attrs.parquet")
artifact.describe()
