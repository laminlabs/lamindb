import lamindb as ln

schema = ln.Schema(
    name="Mini immuno schema",
    features=[
        ln.Feature.get(name="perturbation"),
        ln.Feature.get(name="cell_type_by_model"),
        ln.Feature.get(name="assay_oid"),
        ln.Feature.get(name="donor"),
        ln.Feature.get(name="concentration"),
        ln.Feature.get(name="treatment_time_h"),
    ],
    flexible=True,  # _additional_ columns in a dataframe are validated & annotated
).save()
