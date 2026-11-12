---
execute_via: python
---

# Query tables in storage

This guide walks through querying tabular datasets stored in parquet and related file formats.
Queries stream directly from disk or cloud storage with PyArrow, Polars or DuckDB.

For an in-depth example and a complex dataset, read this blog post: [blog.lamin.ai/1000genomes](https://blog.lamin.ai/1000genomes).

## Experimental: write a LazyFrame without a cache copy

Use `Artifact.from_lazyframe()` to write a Polars LazyFrame directly
to managed S3 storage (or to local storage without an additional cache file).
Pass the frame and a saved schema to validate **inside the sink pipeline**:

```{code-block} python
lazy = ln.Artifact.from_lazyframe(
    lazy_df,
    key="results/counts.parquet",
    schema=my_schema,
)
artifact = lazy.save()
```

Like `from_dataframe()`, this defaults to Parquet and accepts `key`, `description`,
`run`, `revises`, `schema`, `features`, and `parquet_kwargs`. File versions are
non-overwriting by default; no suffix or overwrite flag is needed. CSV output
is not yet supported by this streaming prototype. Polars schema restrictions
match `from_dataframe()`: flexible schemas, including `"valid_features"`, are
not supported. Omitting `key` creates an independent version family unless
`revises` is supplied.

`.save()` evaluates the source once. Each incoming batch is validated with the
existing dataframe curator, and its validated/coerced values are passed onward
to the Parquet sink. Schema planning can read input metadata; validation of
in-memory batches can execute queries on those batches, but does not rescan the
original source. Categorical annotations are combined across batches, subject to
the normal annotation limits.

Choose how dataframe validation failures are handled with `on_validation_error`:

- `"raise"` (default) aborts at the first failing batch, removes partial output,
  and raises the validation error. No artifact is registered, existing version
  heads remain unchanged, and the UID reservation is retained as a tombstone.
  If cleanup fails, the validation error is raised with the cleanup failure as
  its cause; call `lazy.discard()` to retry cleanup.
- `"keep"` completes the file with all original rows, including invalid values.
  If validation fails, the new artifact has `schema=None`, no partial dataframe
  annotations, and a failure report in `artifact.extra_data["validation"]`
  containing `status`, the attempted `schema_uid`, and the first error.
  It is registered as a normal new version and can become the latest version.
  Keep mode rejects `Schema.coerce` or dataframe `Feature.coerce` before any
  source execution or UID reservation: an unconvertible later value cannot be
  preserved in a column whose earlier batches were already coerced.

```{code-block} python
lazy = ln.Artifact.from_lazyframe(
    lazy_df,
    key="results/counts.parquet",
    schema=my_non_coercing_schema,
    on_validation_error="keep",
)
artifact = lazy.save()
```

In keep mode, the first validation failure emits a warning and subsequent
batches can pass through without further validation. `lazy.validation_error`
exposes the first error, bounded to 4096 characters, or `None` if validation has
not failed. Successful validation still attaches the schema and annotations.
Source, I/O, and configuration errors are not treated as validation failures:
they still raise, and an incomplete sink may need `lazy.discard()`.

The existing `Artifact.from_lazy()` flow remains caller-managed. You can keep
writing to `lazy.path` yourself and opt into UID reservation
and same-key deduplication explicitly with `deduplicate="key"`:

```{code-block} python
lazy = ln.Artifact.from_lazy(
    suffix=".parquet",
    overwrite_versions=False,
    key="results/counts.parquet",
    deduplicate="key",
)
lazy_df.sink_parquet(lazy.path.as_posix())
artifact = lazy.save()
```

Without a supplied frame, Lamin cannot validate dataframe content during the
sink; passing a dataframe-content schema in that form is rejected rather than
marking the file validated. Legacy external-metadata schemas remain supported
without a frame.

This opt-in prototype requires the new UID-reservation database migration. Creating
the lazy artifact reserves a normal artifact UID before exposing `lazy.path`,
keeping the existing physical filename
`.lamindb/{uid}.parquet`, virtual key semantics, and `_real_key=None`. No rename or
bucket-side copy is needed. Ordinary artifact version allocation skips reserved
UIDs, including abandoned reservations.
An ordinary unsaved artifact constructed before a conflicting reservation must
be reconstructed to obtain a new UID; its save fails before uploading.

Deduplication checks storage, virtual key, and hash **after** the upload.
Identical content under different keys is allowed. A matching artifact published
while streaming is reused without changing the current version head, and the
redundant object is removed. S3 hashes come from object metadata and may be
multipart ETags, not whole-file MD5 hashes; they are not guaranteed to match a
local-file hash for the same bytes.
Deduplication never overwrites an existing artifact's schema or metadata with
the current attempt's failure report; inspect `lazy.validation_error` for that
attempt.

If a different version is published during streaming, finalization raises a
conflict and retains the uploaded object at `lazy.path`. Call
`lazy.discard()` to remove that object before starting another write. Close the
writer before calling `.save()`. Transient metadata or database failures can be
retried with `.save()` without repeating the upload. UID reservations are
retained as tombstones, so version numbers may have gaps.

This prototype supports the current instance only, single-file Parquet, and
non-overwriting virtual keys. Create the lazy artifact and run `.save()` outside
an enclosing database transaction. It does not provide automatic cleanup or process-restart recovery;
reservations track unfinished paths for manual recovery. S3 deployments should
configure lifecycle cleanup for abandoned multipart uploads. Folder/partitioned
output and other cloud providers are not supported in the opt-in flow. The
existing `from_lazy()` behavior, including folder artifacts, is unchanged when
`deduplicate` is omitted.

```python
import lamindb as ln

db = ln.DB("laminlabs/lamindata")
```

## Open a dataset for queries

To open an **artifact** as a [pyarrow dataset](https://arrow.apache.org/docs/python/dataset.html) call `.open()` with default arguments:

```python
artifact = db.Artifact.filter(
    key__startswith="example_datasets/small", suffix=".parquet"
).first()
dataset = artifact.open()
dataset
```

To open a **collection** (a sharded dataset), call:

```python
collection = db.Collection.get(key="sharded_parquet_collection")
dataset = collection.open()
dataset
```

We recommended to use {class}`~lamindb.Collection` to manage large-scale datasets as this gives you typical lakehouse guarantees. For more background, read {doc}`architecture`.

<!-- #region -->

:::{dropdown} Open sets of artifacts ad hoc

You can also filter artifacts by a prefix and open the resulting `ArtifactSet`:

```python
dataset = db.Artifact.filter(
    key__startswith="example_datasets/small", suffix=".parquet", is_latest=True
).open()
dataset
```

:::

<!-- #endregion -->

## Queries

Filters push down into Parquet row groups — only matching data is read from storage, not the whole file.

```python
import pyarrow.compute as pc

dataset.to_table(filter=pc.field("cell_type").is_valid()).to_pandas().head(5)
```

You can build up from there. Materialize the filtered result once, then compute against it in memory — for example, count rows per cell type without re-reading storage:

```python
counts = (
    dataset.to_table(filter=pc.field("cell_type").is_valid())
    .group_by("cell_type")
    .aggregate([("cell_type", "count")])
    .to_pandas()
)
counts.head()
```

These filters run the same way on every engine below. For heavier patterns — schema evolution, time travel, appends, and concurrent writers — the [LaminDB lakehouse benchmark](https://lamin.ai/blog/lakehouse-benchmark) walks through all six operations end to end across the five engines on a shared genomics dataset.

## Choosing a query engine

By default `Artifact.open()` and `Collection.open()` use `pyarrow` to lazily open dataframes. `polars` can also be used by passing `engine="polars"`. Note also that `.open(engine="polars")` returns a context manager with [LazyFrame](https://docs.pola.rs/api/python/stable/reference/lazyframe/index.html).

```python
with collection.open(engine="polars") as lazy_df:
    display(lazy_df.collect().to_pandas())
```

LaminDB artifacts are plain Parquet files on S3, so any engine that reads Parquet works. DuckDB, Iceberg, and LanceDB are also supported — you pick based on your query pattern, not on how the data was stored.

::::::{tab-set}
:::::{tab-item} PyArrow
`.open()` returns a lazy PyArrow dataset backed by S3. Filters push down into Parquet row groups.

```{code-block} python
import pyarrow.compute as pc

dataset = collection.open()
result = dataset.to_table(filter=pc.field("cell_type").is_valid()).to_pandas()
```

:::::

:::::{tab-item} Polars
`.open(engine="polars")` returns a context manager yielding a Polars LazyFrame backed by S3. No data is read until `.collect()` is called.

```{code-block} python
import polars as pl

with collection.open(engine="polars") as lazy_df:
    result = lazy_df.filter(pl.col("cell_type").is_not_null()).collect()
```

:::::

:::::{tab-item} DuckDB
DuckDB reads Parquet files directly via `read_parquet()`. Use `.cache()` to get a local path, or register a view over S3 paths via `httpfs` for a collection:

```{code-block} python
import duckdb

con = duckdb.connect()
con.execute("INSTALL httpfs; LOAD httpfs;")
con.execute("CREATE OR REPLACE SECRET s3 (TYPE s3, PROVIDER credential_chain);")

s3_paths = [str(a.path) for a in collection.ordered_artifacts.all()]
con.execute(f"CREATE OR REPLACE VIEW data AS SELECT * FROM read_parquet({s3_paths})")

result = con.execute("SELECT * FROM data WHERE cell_type IS NOT NULL").df()
```

:::::

:::::{tab-item} Iceberg
Iceberg requires a one-time ingestion into a catalog-managed table on S3. After ingestion you get native partition pruning, metadata-only schema evolution, and snapshot-based time travel.

```{code-block} python
from pyiceberg.catalog.sql import SqlCatalog
from pyiceberg.expressions import NotNull

arrow = collection.open().to_table()

catalog = SqlCatalog(
    "local", uri="sqlite:///iceberg_catalog.db", warehouse="/tmp/iceberg_warehouse"
)
if not catalog.namespace_exists("demo"):
    catalog.create_namespace("demo")
table = catalog.create_table("demo.data", schema=arrow.schema)
table.append(arrow)

result = table.scan(row_filter=NotNull("cell_type")).to_arrow().to_pandas()
```

:::::

:::::{tab-item} LanceDB
LanceDB ingests data into Lance columnar format on S3 — the only engine here that copies data out of the source Parquet files. In exchange you get combined SQL filtering and vector search, plus versioned appends with time travel.

```{code-block} python
import lancedb
import pyarrow as pa

arrow = collection.open().to_table()
schema = pa.schema([
    f.with_type(pa.string()) if pa.types.is_dictionary(f.type) else f
    for f in arrow.schema
])
arrow = arrow.cast(schema)
db_lance = lancedb.connect("/tmp/lancedb_warehouse")
table = db_lance.create_table("data", data=arrow, mode="overwrite")

result = table.search().where("cell_type IS NOT NULL", prefilter=True).to_pandas()
```

:::::
::::::
