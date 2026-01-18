import lamindb as ln
from datetime import date

df = ln.examples.datasets.mini_immuno.get_dataset1(otype="DataFrame")

temperature = ln.Feature(name="temperature", dtype=float).save()
date_of_study = ln.Feature(name="date_of_study", dtype=date).save()
external_schema = ln.Schema(features=[temperature, date_of_study]).save()

concentration = ln.Feature(name="concentration", dtype=str).save()
donor = ln.Feature(name="donor", dtype=str, nullable=True).save()
schema = ln.Schema(
    features=[concentration, donor],
    slots={"__external__": external_schema},
    otype="DataFrame",
).save()

artifact = ln.Artifact.from_dataframe(
    df,
    key="examples/dataset1.parquet",
    features={"temperature": 21.6, "date_of_study": date(2024, 10, 1)},
    schema=schema,
).save()
artifact.describe()
