import lamindb as ln

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
    spatialdata, key="examples/spatialdata1.zarr", schema=sdata_schema
).save()
artifact.describe()
