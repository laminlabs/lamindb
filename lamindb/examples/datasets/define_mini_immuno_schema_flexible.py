import lamindb as ln

schema = ln.Schema(
    name="mini_immuno",
    features=[  # required features
        ln.Feature.get(name="perturbation"),
        ln.Feature.get(name="cell_type_by_model"),
        ln.Feature.get(name="assay_oid"),
    ],
    itype=ln.Feature.get(name="ImmunoSamples"),
    flexible=True,  # additional columns in a dataframe are validated & annotated if they match the feature type
).save()
