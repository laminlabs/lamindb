import lamindb as ln
import bionty as bt


attrs_schema = ln.Schema(
    features=[
        ln.Feature(name="bio", dtype=dict).save(),
        ln.Feature(name="tech", dtype=dict).save(),
    ],
).save()

sample_schema = ln.Schema(
    features=[
        ln.Feature(name="disease", dtype=bt.Disease, coerce_dtype=True).save(),
        ln.Feature(
            name="developmental_stage",
            dtype=bt.DevelopmentalStage,
            coerce_dtype=True,
        ).save(),
    ],
).save()

tech_schema = ln.Schema(
    features=[
        ln.Feature(name="assay", dtype=bt.ExperimentalFactor, coerce_dtype=True).save(),
    ],
).save()

obs_schema = ln.Schema(
    features=[
        ln.Feature(name="sample_region", dtype="str").save(),
    ],
).save()

# Schema enforces only registered Ensembl Gene IDs are valid (maximal_set=True)
varT_schema = ln.Schema(itype=bt.Gene.ensembl_gene_id, maximal_set=True).save()

sdata_schema = ln.Schema(
    name="spatialdata_blobs_schema",
    otype="SpatialData",
    slots={
        "attrs:bio": sample_schema,
        "attrs:tech": tech_schema,
        "attrs": attrs_schema,
        "tables:table:obs": obs_schema,
        "tables:table:var.T": varT_schema,
    },
).save()
