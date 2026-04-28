---
execute_via: python
---

# Stream datasets from storage

This guide walks through streaming datasets from disk or cloud storage.

```python
# replace with your username and S3 bucket
!lamin login testuser1
!lamin init --storage s3://lamindb-ci/test-arrays
```

Import lamindb and track this notebook.

```python
import lamindb as ln
import numpy as np

ln.track()
db = ln.Artifact.connect("laminlabs/lamindata")  # we'll pull dataset from there
```

## DataFrame

### Streaming from a single artifact

A dataframe stored as sharded `parquet`.

```python
artifact = db.Artifact.get(key="sharded_parquet")
```

```python
artifact.path.view_tree()
```

```python
dataset = artifact.open()
```

This returns a [pyarrow dataset](https://arrow.apache.org/docs/python/dataset.html).

```python
dataset
```

```python
dataset.head(5).to_pandas()
```

### Streaming from a set of artifacts

You can open several parquet files as a single dataset by calling `.open()` on the result of a query:

```python
dataset = db.Artifact.filter(key__startswith="example_datasets/small", suffix=".parquet").open()  # open an ArtifactSet for streaming
dataset
```

The same is possible for the artifacts in a collection:

```python
collection = db.Collection.get(key="sharded_parquet_collection")
dataset = collection.open()
dataset
```

Once you have a storage-backed dataset, you can query it like this:

```python
dataset.to_table().to_pandas()
```

By default `Artifact.open()` and `Collection.open()` use `pyarrow` to lazily open dataframes. `polars` can be also used by passing `engine="polars"`. Note also that `.open(engine="polars")` returns a context manager with [LazyFrame](https://docs.pola.rs/api/python/stable/reference/lazyframe/index.html).

```python
with collection.open(engine="polars", use_fsspec=True) as lazy_df:
    display(lazy_df.collect().to_pandas())
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

This object is an `AnnDataAccessor` object, an `AnnData` object backed in the cloud:

```python
adata
```

Without subsetting, the `AnnDataAccessor` object references underlying lazy `h5` or `zarr` arrays:

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

It is also possible to add columns to `.obs` and `.var` of cloud AnnData objects without downloading them. First, create a new `AnnData` `zarr` artifact:

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

It is also possible to access `AnnData` objects inside `SpatialData` `tables`:

```python
artifact = ln.Artifact.connect("laminlabs/lamindata").get(
    key="visium_aligned_guide_min.zarr"
)

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

<!-- #region -->

```python
adata = table_subset.to_memory()
```

<!-- #endregion -->

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
