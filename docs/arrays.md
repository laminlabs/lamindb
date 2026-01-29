# Slice & stream arrays

We saw how LaminDB allows to query & search across artifacts using registries: {doc}`registries`.

Let us now query the datasets in storage themselves. Here, we show how to subset an `AnnData` and generic `HDF5` and `zarr` collections accessed in the cloud.

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
```

We'll need some test data:

```python
ln.Artifact("s3://lamindb-ci/test-arrays/pbmc68k.h5ad").save()
ln.Artifact("s3://lamindb-ci/test-arrays/testfile.hdf5").save()
```

## AnnData

An `h5ad` artifact stored on s3:

```python
artifact = ln.Artifact.get(key="pbmc68k.h5ad")
```

```python
artifact.path
```

```python
access = artifact.open()
```

This object is an `AnnDataAccessor` object, an `AnnData` object backed in the cloud:

```python
access
```

Without subsetting, the `AnnDataAccessor` object references underlying lazy `h5` or `zarr` arrays:

```python
access.X
```

You can subset it like a normal `AnnData` object:

```python
obs_idx = access.obs.cell_type.isin(["Dendritic cells", "CD14+ Monocytes"]) & (
    access.obs.percent_mito <= 0.05
)
access_subset = access[obs_idx]
access_subset
```

Subsets load arrays into memory upon direct access:

```python
access_subset.X
```

To load the entire subset into memory as an actual `AnnData` object, use `to_memory()`:

```python
adata_subset = access_subset.to_memory()

adata_subset
```

### Add a column to a cloud AnnData object

It is also possible to add columns to `.obs` and `.var` of cloud AnnData objects without downloading them.

Create a new `AnnData` `zarr` artifact.

```python
adata_subset.write_zarr("adata_subset.zarr")
```

```python
artifact = ln.Artifact(
    "adata_subset.zarr", description="test add column to adata"
).save()
```

```python
artifact
```

```python
with artifact.open(mode="r+") as access:
    access.add_column(where="obs", col_name="ones", col=np.ones(access.shape[0]))
    display(access)
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

## Parquet

A dataframe stored as sharded `parquet`.

Note that it is also possible to register and access Hugging Face paths. For this `huggingface_hub` package should be installed.

```python
artifact = ln.Artifact.connect("laminlabs/lamindata").get(key="sharded_parquet")
```

```python
artifact.path.view_tree()
```

```python
backed = artifact.open()
```

This returns a [pyarrow dataset](https://arrow.apache.org/docs/python/dataset.html).

```python
backed
```

```python
backed.head(5).to_pandas()
```

It is also possible to open a collection of cloud artifacts.

```python
collection = ln.Collection.connect("laminlabs/lamindata").get(
    key="sharded_parquet_collection"
)
```

```python
backed = collection.open()
```

```python
backed
```

```python
backed.to_table().to_pandas()
```

By default `Artifact.open()` and `Collection.open()` use `pyarrow` to lazily open dataframes. `polars` can be also used by passing `engine="polars"`. Note also that `.open(engine="polars")` returns a context manager with [LazyFrame](https://docs.pola.rs/api/python/stable/reference/lazyframe/index.html).

```python
with collection.open(engine="polars", use_fsspec=True) as lazy_df:
    display(lazy_df.collect().to_pandas())
```

Yet another way to open several parquet files as a single dataset is via calling `.open()` directly for a query set.

```python
backed = ln.Artifact.filter(suffix=".parquet").open()
```

```python
backed
```

```python
# clean up test instance
ln.setup.delete("test-arrays", force=True)
```
