---
execute_via: python
---

# Query tables in storage

This guide walks through querying tabular datasets stored in parquet and related file formats.
Queries stream directly from disk or cloud storage with PyArrow, Polars or DuckDB.

```python
import lamindb as ln

db = ln.DB("laminlabs/lamindata")
```

## Stream a table from storage

Start with a single parquet file. Get the artifact and open it — nothing is downloaded, you get a lazy [pyarrow dataset](https://arrow.apache.org/docs/python/dataset.html):

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

The same `.open()` call works across parquet files. Call it on a query to stream a whole set of Parquet artifacts as one dataset:

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
