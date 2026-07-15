---
execute_via: python
---

# Stream arrays from storage

This guide covers streaming array-like data — AnnData, SpatialData, and generic HDF5 — directly from disk or cloud storage. For tabular data (DataFrame, Parquet, and querying with PyArrow, Polars, DuckDB, Iceberg, or LanceDB), see {doc}`Query tables <parquet>`.

```bash
# replace with your username and S3 bucket
lamin login testuser1
lamin init --storage s3://lamindb-ci/test-arrays
```

Import lamindb and track this notebook.

```python
import lamindb as ln
import numpy as np

db = ln.DB("laminlabs/lamindata")  # we'll pull the SpatialData example from there
```

## AnnData

We'll need some test data:

```python
ln.Artifact("s3://lamindb-ci/test-arrays/pbmc68k.h5ad").save()
ln.Artifact("s3://lamindb-ci/test-arrays/testfile.hdf5").save()
```

An `h5ad` artifact stored on s3:

```python
artifact = ln.Artifact.get(key="pbmc68k.h5ad")
```

```python
artifact.path
```

```python
adata = artifact.open()
```

This is an `AnnDataAccessor` — an `AnnData` backed by cloud storage:

```python
adata
```

Without subsetting, it references the underlying lazy `h5` or `zarr` arrays:

```python
adata.X
```

You can subset it like a normal `AnnData` object:

```python
obs_idx = adata.obs.cell_type.isin(["Dendritic cells", "CD14+ Monocytes"]) & (
    adata.obs.percent_mito <= 0.05
)
adata_subset = adata[obs_idx]
adata_subset
```

Subsets load arrays into memory upon direct access:

```python
adata_subset.X
```

To load the entire subset into memory as an actual `AnnData` object, use `to_memory()`:

```python
adata_subset.to_memory()
```

You can also add columns to `.obs` and `.var` of a cloud `AnnData` without downloading it. First, create a new `AnnData` `zarr` artifact:

```python
adata_subset.to_memory().write_zarr("adata_subset.zarr")
artifact = ln.Artifact(
    "adata_subset.zarr", description="test add column to adata"
).save()
```

This is how you add a column:

```python
with artifact.open(mode="r+") as adata_accessor:
    adata_accessor.add_column(where="obs", col_name="ones", col=np.ones(adata_accessor.shape[0]))
    display(adata_accessor)
```

The version of the artifact is updated after the modification.

```python
artifact
```

```python
artifact.delete(permanent=True)
```

## SpatialData

You can also access `AnnData` objects inside `SpatialData` tables:

```python
artifact = db.Artifact.get(key="visium_aligned_guide_min.zarr")

access = artifact.open()
```

```python
access
```

```python
access.tables
```

This gives you the same `AnnDataAccessor` object as for a normal `AnnData`.

```python
table = access.tables["table"]

table
```

You can subset it and read into memory as an actual `AnnData`:

```python
table_subset = table[table.obs["clone"] == "diploid"]

table_subset
```

```python
adata = table_subset.to_memory()
```

## Generic HDF5

Let us query a generic HDF5 artifact:

```python
artifact = ln.Artifact.get(key="testfile.hdf5")
```

And get a backed accessor:

```python
backed = artifact.open()
```

The returned object contains the `.connection` and `h5py.File` or `zarr.Group` in `.storage`

```python
backed
```

```python
backed.storage
```

```python
# clean up test instance
ln.setup.delete("test-arrays", force=True)
```
