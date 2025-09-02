import lamindb as ln

study_metadata_schema = ln.Schema(
    features=[
        ln.Feature(name="temperature", dtype=float).save(),
        ln.Feature(name="experiment_id", dtype=str).save(),
    ],
).save()
