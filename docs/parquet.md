---
execute_via: python
---

# Query tables in storage

This guide walks through querying tabular data stored as Parquet — streaming directly from disk or cloud storage with PyArrow, Polars, DuckDB, Iceberg, or LanceDB.

Import lamindb and track this notebook.

```python
import lamindb as ln

ln.track()  # this step is optional
db = ln.DB("laminlabs/lamindata")  # we'll pull datasets from there
```

## Stream a table from storage

Start with a single Parquet file. Get the artifact and open it — nothing is downloaded, you get a lazy [pyarrow dataset](https://arrow.apache.org/docs/python/dataset.html) backed by storage.

```python
artifact = db.Artifact.filter(
    key__startswith="example_datasets/small", suffix=".parquet"
).first()

dataset = artifact.open()
dataset
```

Peek at the first few rows:

```python
dataset.head(5).to_pandas()
```

The same `.open()` call scales to many files. Call it on a query to stream a whole set of Parquet artifacts as one dataset:

```python
dataset = db.Artifact.filter(
    key__startswith="example_datasets/small", suffix=".parquet", is_latest=True
).open()
dataset
```

Or on a collection, which streams all of its member Parquet files together:

```python
collection = db.Collection.get(key="sharded_parquet_collection")
dataset = collection.open()
dataset.to_table().to_pandas()
```

Whether the data is one file or a thousand shards, you query it the same way.

## Queries

Filters push down into Parquet row groups — only matching data is read from storage, not the whole file.

```python
import pyarrow.compute as pc

dataset.to_table(filter=pc.field("disease").is_valid()).to_pandas().head(5)
```

You can build up from there. Materialize the filtered result once, then compute against it in memory — for example, count rows per disease without re-reading storage:

```python
counts = (
    dataset.to_table(filter=pc.field("disease").is_valid())
    .group_by("disease")
    .aggregate([("disease", "count")])
    .to_pandas()
)
counts.head()
```

These filters run the same way on every engine below. For heavier patterns — schema evolution, time travel, appends, and concurrent writers — the [LaminDB lakehouse benchmark](https://lamin.ai/blog/lakehouse-benchmark) walks through all six operations end to end across the five engines on a shared genomics dataset.

## Choosing a query engine

By default `Artifact.open()` and `Collection.open()` use `pyarrow` to lazily open dataframes. `polars` can also be used by passing `engine="polars"`. Note also that `.open(engine="polars")` returns a context manager with [LazyFrame](https://docs.pola.rs/api/python/stable/reference/lazyframe/index.html).

```python
with collection.open(engine="polars", use_fsspec=True) as lazy_df:
    display(lazy_df.collect().to_pandas())
```

LaminDB artifacts are plain Parquet files on S3, so any engine that reads Parquet works. DuckDB, Iceberg, and LanceDB are also supported — you pick based on your query pattern, not on how the data was stored.

::::::{tab-set}
:::::{tab-item} PyArrow
`.open()` returns a lazy PyArrow dataset backed by S3. Filters push down into Parquet row groups.

```python
import pyarrow.compute as pc

dataset = collection.open()
result = dataset.to_table(filter=pc.field("disease").is_valid()).to_pandas()
```

:::::

:::::{tab-item} Polars
`.open(engine="polars")` returns a context manager yielding a Polars LazyFrame backed by S3. No data is read until `.collect()` is called.

```python
import polars as pl

with collection.open(engine="polars") as lazy_df:
    result = lazy_df.filter(pl.col("disease").is_not_null()).collect()
```

:::::

:::::{tab-item} DuckDB
DuckDB reads Parquet files directly via `read_parquet()`. Use `.cache()` to get a local path, or register a view over S3 paths via `httpfs` for a collection:

```python
import duckdb

con = duckdb.connect()
con.execute("INSTALL httpfs; LOAD httpfs;")
con.execute("CREATE OR REPLACE SECRET s3 (TYPE s3, PROVIDER credential_chain);")

s3_paths = [str(a.path) for a in collection.ordered_artifacts.all()]
con.execute(f"CREATE OR REPLACE VIEW data AS SELECT * FROM read_parquet({s3_paths})")

result = con.execute("SELECT * FROM data WHERE disease IS NOT NULL").df()
```

:::::

:::::{tab-item} Iceberg
Iceberg requires a one-time ingestion into a catalog-managed table on S3. After ingestion you get native partition pruning, metadata-only schema evolution, and snapshot-based time travel.

```python
from pyiceberg.catalog.sql import SqlCatalog
from pyiceberg.expressions import NotNull

arrow = collection.open().to_table()

catalog = SqlCatalog(
    "local", uri="sqlite:///iceberg_catalog.db", warehouse="s3://your-bucket/iceberg_warehouse"
)
catalog.create_namespace("demo")
table = catalog.create_table("demo.data", schema=arrow.schema)
table.append(arrow)

result = table.scan(row_filter=NotNull("disease")).to_arrow().to_pandas()
```

:::::

:::::{tab-item} LanceDB
LanceDB ingests data into Lance columnar format on S3 — the only engine here that copies data out of the source Parquet files. In exchange you get combined SQL filtering and vector search, plus versioned appends with time travel.

```python
import lancedb

arrow = collection.open().to_table()
db_lance = lancedb.connect("s3://your-bucket/lancedb_warehouse")
table = db_lance.create_table("data", data=arrow, mode="overwrite")

result = table.search().where("disease IS NOT NULL", prefilter=True).to_pandas()
```

:::::
::::::

## Performance

The charts below compare PyArrow, Polars, DuckDB, Iceberg, and LanceDB on the same LaminDB collection (8,929 rows, 6 Parquet shards) on S3. See the [full blog post](https://lamin.ai/blog/lakehouse-benchmark) for methodology and caveats.

<!-- PLOT: setup_cost.svg -->

![Setup cost](https://lamin-site-assets.s3.amazonaws.com/.lamindb/Lf8f0LJY63quZ3n70000.svg)

<!-- PLOT: query_times.svg -->

![Query times](https://lamin-site-assets.s3.amazonaws.com/.lamindb/d2r3p1yUGrcVTLtw0000.svg)

<!-- PLOT: write_path.svg -->

![Write-path times](https://lamin-site-assets.s3.amazonaws.com/.lamindb/VnVruqKX9KK0uhUw0000.svg)

```python
# clean up test instance
ln.setup.delete("test-tables", force=True)
```
