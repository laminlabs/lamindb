---
execute_via: python
---

# Stream datasets from storage

This guide walks through streaming datasets from disk or cloud storage.

```
# replace with your username and S3 bucket
!lamin login testuser1
!lamin init --storage s3://lamindb-ci/test-arrays
```

Import lamindb and track this notebook.

```
import lamindb as ln
import numpy as np

ln.track()
db = ln.DB("laminlabs/lamindata")  # we'll pull dataset from there
```

For tabular data — DataFrame, Parquet, querying with PyArrow, Polars, DuckDB, Iceberg, or LanceDB — see [Query tables](tables.md). This guide covers array-like data: AnnData, SpatialData, and generic HDF5.

## AnnData

We'll need some test data:

```
ln.Artifact("s3://lamindb-ci/test-arrays/pbmc68k.h5ad").save()
ln.Artifact("s3://lamindb-ci/test-arrays/testfile.hdf5").save()
```

An `h5ad` artifact stored on s3:

```
artifact = ln.Artifact.get(key="pbmc68k.h5ad")
```

```
artifact.path
```

```
adata = artifact.open()
```

This object is an `AnnDataAccessor` object, an `AnnData` object backed in the cloud:

```
adata
```

Without subsetting, the `AnnDataAccessor` object references underlying lazy `h5` or `zarr` arrays:

```
adata.X
```

You can subset it like a normal `AnnData` object:

```
obs_idx = adata.obs.cell_type.isin(["Dendritic cells", "CD14+ Monocytes"]) & (
    adata.obs.percent_mito <= 0.05
)
adata_subset = adata[obs_idx]
adata_subset
```

Subsets load arrays into memory upon direct access:

```
adata_subset.X
```

To load the entire subset into memory as an actual `AnnData` object, use `to_memory()`:

```
adata_subset.to_memory()
```

It is also possible to add columns to `.obs` and `.var` of cloud AnnData objects without downloading them. First, create a new `AnnData` `zarr` artifact:

```
adata_subset.to_memory().write_zarr("adata_subset.zarr")
artifact = ln.Artifact(
    "adata_subset.zarr", description="test add column to adata"
).save()
```

This is how you add a column:

```
with artifact.open(mode="r+") as adata_accessor:
    adata_accessor.add_column(where="obs", col_name="ones", col=np.ones(adata_accessor.shape[0]))
    display(adata_accessor)
```

The version of the artifact is updated after the modification.

```
artifact
```

```
artifact.delete(permanent=True)
```

## SpatialData

It is also possible to access `AnnData` objects inside `SpatialData` `tables`:

```
artifact = ln.Artifact.connect("laminlabs/lamindata").get(
    key="visium_aligned_guide_min.zarr"
)

access = artifact.open()
```

```
access
```

```
access.tables
```

This gives you the same `AnnDataAccessor` object as for a normal `AnnData`.

```
table = access.tables["table"]

table
```

You can subset it and read into memory as an actual `AnnData`:

```
table_subset = table[table.obs["clone"] == "diploid"]

table_subset
```

```
adata = table_subset.to_memory()
```

## Generic HDF5

Let us query a generic HDF5 artifact:

```
artifact = ln.Artifact.get(key="testfile.hdf5")
```

And get a backed accessor:

```
backed = artifact.open()
```

The returned object contains the `.connection` and `h5py.File` or `zarr.Group` in `.storage`

```
backed
```

```
backed.storage
```

```
# clean up test instance
ln.setup.delete("test-arrays", force=True)
```