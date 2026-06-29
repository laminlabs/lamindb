---
execute_via: python
---

# Query tables

This guide walks through querying tabular Parquet data from LaminDB storage.

```python
# replace with your username and S3 bucket
!lamin login testuser1
!lamin init --storage s3://your-bucket/test-tables
```

```python
import lamindb as ln
import pandas as pd
import pyarrow.compute as pc

ln.track()
```

Save a DataFrame as a versioned, lineage-tracked artifact:

```python
df = pd.DataFrame({
    "chrom": ["chr1", "chr17", "chr17"],
    "pos":   [100_000, 7_600_000, 7_700_000],
    "gene":  ["BRCA2", "TP53", "BRCA1"],
})
artifact = ln.Artifact.from_dataframe(df, key="demo/variants.parquet").save()
```

Query it back — streams from S3, filter pushes down into Parquet row groups:

```python
result = (
    artifact.open()
    .to_table(filter=pc.field("chrom") == "chr17")
    .to_pandas()
)
result
```

## PyArrow

`artifact.open()` returns a lazy [PyArrow dataset](https://arrow.apache.org/docs/python/dataset.html) backed by S3. Filters push down into Parquet row groups — only matching data is read.

```python
artifact = ln.Artifact.get(key="demo/variants.parquet")
dataset = artifact.open()

result = dataset.to_table(
    filter=(pc.field("chrom") == "chr17"),
    columns=["chrom", "pos", "gene"],
).to_pandas()
result
```

For a **collection** of many Parquet shards (e.g. one file per sample), call `.open()` on the collection:

```python
collection = ln.Collection.get(key="sharded_parquet_collection")
dataset = collection.open()   # opens all shards as one unified dataset
dataset.to_table(filter=pc.field("chrom") == "chr17").to_pandas()
```

## Polars

`artifact.open(engine="polars")` returns a context manager yielding a Polars [LazyFrame](https://docs.pola.rs/api/python/stable/reference/lazyframe/index.html) backed by S3. No data is read until `.collect()` is called.

```python
import polars as pl

artifact = ln.Artifact.get(key="demo/variants.parquet")

with artifact.open(engine="polars") as lazy_df:
    result = (
        lazy_df
        .filter(pl.col("chrom") == "chr17")
        .select(["chrom", "pos", "gene"])
        .collect()
    )
result
```

## DuckDB

DuckDB registers a view over S3 paths via `httpfs`. The view is lazy — nothing is read until a SQL query is issued.

```python
import duckdb

con = duckdb.connect()
con.execute("INSTALL httpfs; LOAD httpfs;")
con.execute("CREATE OR REPLACE SECRET s3 (TYPE s3, PROVIDER credential_chain);")

collection = ln.Collection.get(key="sharded_parquet_collection")
s3_paths = [str(a.path) for a in collection.ordered_artifacts.all()]
con.execute(f"CREATE OR REPLACE VIEW variants AS SELECT * FROM read_parquet({s3_paths})")

con.execute("SELECT chrom, pos, gene FROM variants WHERE chrom = 'chr17'").df()
```

For a single artifact, use `.cache()` to get a local path:

```python
path = artifact.cache()
duckdb.sql(f"SELECT * FROM read_parquet('{path}') WHERE chrom = 'chr17'").df()
```

## Iceberg

Iceberg requires a one-time ingestion into a catalog-managed table. After ingestion, queries use native partition pruning, schema evolution is metadata-only, and any historical snapshot is addressable by ID.

```python
from pyiceberg.catalog.sql import SqlCatalog
from pyiceberg.expressions import EqualTo
from pyiceberg.types import BooleanType

arrow = artifact.open().to_table()

catalog = SqlCatalog(
    "local",
    uri="sqlite:///iceberg_catalog.db",
    warehouse="s3://your-bucket/iceberg_warehouse",
)
catalog.create_namespace("genomics")
table = catalog.create_table("genomics.variants", schema=arrow.schema)
table.append(arrow)
```

Query with native pushdown:

```python
table.scan(row_filter=EqualTo("chrom", "chr17")).to_arrow().to_pandas()
```

Schema evolution — add a column with no Parquet rewrite:

```python
with table.update_schema() as update:
    update.add_column("qc_pass", BooleanType())
```

Time travel — read any historical snapshot:

```python
first_snapshot = table.history()[0].snapshot_id
table.scan(snapshot_id=first_snapshot).to_arrow()
```

Note: a SQLite catalog is used here for portability. Production deployments would use a Glue or REST catalog.

## LanceDB

LanceDB ingests data into Lance columnar format, enabling combined SQL filtering and vector search with versioned appends.

```python
import lancedb

arrow = artifact.open().to_table()
db = lancedb.connect("s3://your-bucket/lancedb_warehouse")
table = db.create_table("variants", data=arrow, mode="overwrite")
```

Query with native pushdown:

```python
table.search().where("chrom = 'chr17'", prefilter=True).to_pandas()
```

Append and time travel:

```python
table.add(new_sample_arrow)
table.checkout(1)       # version 1 = pre-append state
table.checkout_latest()
```

## Performance

The charts below are from the [LaminDB lakehouse benchmark](https://lamin.ai/laminlabs/lakehouse-benchmarks) — all engines on the same 1000 Genomes CNV dataset (8,929 rows, 6 Parquet shards) on S3. See the [full blog post](https://lamin.ai/blog/lakehouse-benchmark) for methodology and caveats.

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