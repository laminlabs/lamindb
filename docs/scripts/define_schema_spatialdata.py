import lamindb as ln
import bionty as bt

# a very comprehensive schema for different slots of a SpatialData object

# define or query features
bio_dict = ln.Feature(name="bio", dtype=dict).save()
tech_dict = ln.Feature(name="tech", dtype=dict).save()
disease = ln.Feature(name="disease", dtype=bt.Disease, coerce=True).save()
developmental_stage = ln.Feature(
    name="developmental_stage",
    dtype=bt.DevelopmentalStage,
    coerce=True,
).save()
assay = ln.Feature(name="assay", dtype=bt.ExperimentalFactor, coerce=True).save()
sample_region = ln.Feature(name="sample_region", dtype=str).save()
analysis = ln.Feature(name="analysis", dtype=str).save()

# define or query schema components
attrs_schema = ln.Schema([bio_dict, tech_dict]).save()
sample_schema = ln.Schema([disease, developmental_stage]).save()
tech_schema = ln.Schema([assay]).save()
obs_schema = ln.Schema([sample_region]).save()
uns_schema = ln.Schema([analysis]).save()
# enforces only registered Ensembl Gene IDs pass validation (maximal_set=True)
varT_schema = ln.Schema(itype=bt.Gene.ensembl_gene_id, maximal_set=True).save()

# compose the SpatialData schema
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
