import lamindb as ln

# curate a SpatialData
spatialdata = ln.core.datasets.spatialdata_blobs()
sdata_schema = ln.Schema.get(name="spatialdata_blobs_schema")
curator = ln.curators.SpatialDataCurator(spatialdata, sdata_schema)
try:
    curator.validate()
except ln.errors.ValidationError:
    pass

spatialdata.tables["table"].var.drop(index="ENSG00000999999", inplace=True)

# validate again (must pass now) and save artifact
artifact = ln.Artifact.from_spatialdata(
    spatialdata, key="example_datasets/spatialdata1.zarr", schema=sdata_schema
)
assert artifact.schema == sdata_schema
