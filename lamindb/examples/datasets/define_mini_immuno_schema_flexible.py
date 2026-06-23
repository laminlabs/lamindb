import lamindb as ln

mitype = ln.Feature.get(name="mini_immuno", is_type=True)

schema = ln.Schema(
    name="mini_immuno",
    features=[
        ln.Feature.get(name="perturbation", type=mitype),
        ln.Feature.get(name="cell_type_by_model", type=mitype),
        ln.Feature.get(name="assay_oid", type=mitype),
        ln.Feature.get(name="donor", type=mitype),
        ln.Feature.get(name="concentration", type=mitype),
        ln.Feature.get(name="treatment_time_h", type=mitype),
    ],
    flexible=True,  # _additional_ columns in a dataframe are validated & annotated
).save()
