import lamindb as ln
import bionty as bt


# define sample schema
sample_schema = ln.Schema(
    name="blobs_sample_level_metadata",
    features=[
        ln.Feature(name="assay", dtype=bt.ExperimentalFactor).save(),
        ln.Feature(name="disease", dtype=bt.Disease).save(),
    ],
).save()

# define table obs schema
blobs_obs_schema = ln.Schema(
    name="blobs_obs_level_metadata",
    features=[
        ln.Feature(name="sample_region", dtype="str").save(),
    ],
).save()

# define table var schema
blobs_var_schema = ln.Schema(
    name="blobs_var_schema", itype=bt.Gene.ensembl_gene_id, dtype=int
).save()

# define composite schema
spatialdata_schema = ln.Schema(
    name="blobs_spatialdata_schema",
    otype="SpatialData",
    components={
        "sample": sample_schema,
        "table:obs": blobs_obs_schema,
        "table:var": blobs_var_schema,
    },
).save()

# curate a SpatialData
spatialdata = ln.core.datasets.spatialdata_blobs()
curator = ln.curators.SpatialDataCurator(spatialdata, spatialdata_schema)
try:
    curator.validate()
except ln.errors.ValidationError:
    pass

spatialdata.tables["table"].var.drop(index="ENSG00000999999", inplace=True)

# validate again (must pass now) and save artifact
artifact = ln.Artifact.from_spatialdata(
    spatialdata, key="example_datasets/spatialdata1.zarr", schema=spatialdata_schema
)
assert artifact.schema == spatialdata_schema
