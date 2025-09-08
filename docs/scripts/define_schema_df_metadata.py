import lamindb as ln

study_metadata_schema = ln.Schema(
    name="Study metadata schema",
    features=[
        ln.Feature(name="temperature", dtype=float).save(),
        ln.Feature(name="experiment", dtype=str).save(),
    ],
).save()
