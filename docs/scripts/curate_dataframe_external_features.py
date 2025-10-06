import lamindb as ln

df = ln.examples.datasets.small_dataset1(otype="DataFrame")

species = ln.Feature(name="species", dtype="str").save()
split = ln.Feature(name="split", dtype="str").save()
external_schema = ln.Schema(features=[species, split]).save()

feat1 = ln.Feature(name="feat1", dtype="int").save()
feat2 = ln.Feature(name="feat2", dtype="int").save()
schema = ln.Schema(
    features=[feat1, feat2], slots={"__external__": external_schema}, otype="DataFrame"
).save()

artifact = ln.Artifact.from_dataframe(
    df,
    features={"species": "bird", "split": "train"},
    schema=schema,
    description="test dataframe with external features",
).save()
artifact.describe()
