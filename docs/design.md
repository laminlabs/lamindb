# System design

```{toctree}
:maxdepth: 1
:hidden:

acid
idempotency
```

LaminDB is a distributed data management system like git that can be run or hosted anywhere. It just needs a SQLite or Postgres database and at least one storage location (file system, S3, GCP, Hugging Face, ...). Creating a local LaminDB instance after `pip install lamindb` is as easy as:

::::{tab-set}
:::{tab-item} Shell

```bash
lamin init
```

:::
:::{tab-item} Py
:sync: python

```python
import lamindb as ln
ln.setup.init()
```

:::
:::{tab-item} R
:sync: r

```R
library(laminr)
lamin_init()
```

:::
::::

Or you connect to an existing remote database:

::::{tab-set}
:::{tab-item} Shell

```bash
lamin connect --here account/instance  # --here localizes your connection to the current directory
```

:::
:::{tab-item} Py
:sync: python

```python
import lamindb as ln
ln.connect("account/instance")
```

:::
:::{tab-item} R
:sync: r

```R
library(laminr)
ln <- import_module("lamindb")
ln <- ln$connect("account/instance")
```

:::
::::

For more configuration, see {doc}`docs:setup`. LaminDB instances work standalone but can optionally be managed by LaminHub.

## Lakehouse architecture

Working with a high number of raw files across different sources almost inevitably leads to fragile data organization. This brittleness is amplified when working with agents: they prioritize solving the immediate task over long-term maintainability, they make frequent mistakes, and their concurrent read/write patterns can quickly corrupt a purely file-based architecture. Lakehouse frameworks solve these problems with [ACID transactions](https://en.wikipedia.org/wiki/ACID) to prevent partial writes, with schema enforcement to prevent inconsistent datasets, and with time travel to easily restore erroneous written datasets.
And, as discussed earlier, they also make agents more efficient. So, let's briefly review available options.

### Frameworks

<figure style="float: right; width: 400px; margin-left: 0.5rem">
  <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/OgVhDACCMhzGKC4t0001.svg" />
  *File layout of an Iceberg table.*
</figure>

Today's most popular framework is **Iceberg**.[^apache-iceberg] Like Delta Lake[^delta][^databricks] and Apache Hudi,[^hudi] Iceberg provides ACID transactions and "time travel" by organizing parquet files into snapshots, managed by manifest and metadata files (**Figure 4**). However, this file-based metadata introduces costs (snapshot creation is expensive dictating large, infrequent writes), optimistic concurrency leads to conflicts between simultaneous writers, and coordinating updates on S3 requires an external catalog like AWS Glue or Nessie.[^nessie]

<div style="float: right; width: 65%; margin: 0.5rem 0 1rem 1.5rem; font-size: 0.85em;">

| Feature                                  | Raw S3 | Iceberg | DuckLake | LaminDB |
| ---------------------------------------- | ------ | ------- | -------- | ------- |
| Data lake (file management & annotation) | ✅     | ❌      | ❌       | ✅      |
| ACID transactions                        | ❌     | ✅      | ✅       | ✅ ¹    |
| Time travel / snapshot version isolation | ❌     | ✅      | ✅       | ✅ ²    |
| Schema evolution without rewriting data  | ❌     | ✅ ³    | ✅ ³     | ✅ ³    |
| Write-Audit-Publish workflow             | ❌     | ✅      | ❌       | ✅ ⁴    |
| Query engine independence                | ✅     | ✅      | ❌       | ✅      |
| Concurrent writers                       | ❌ ⁵   | ❌      | ✅       | ✅      |
| Automatic maintenance                    | ❌     | ❌      | ✅ ⁶     | ✅ ⁶    |
| Native multi-table transactions          | ❌     | ❌      | ✅       | ❌      |
| Dataset formats beyond tables            | ✅     | ❌      | ❌       | ✅      |
| Data lineage                             | ❌     | ❌      | ❌       | ✅      |
| Registries/ontologies                    | ❌     | ❌      | ❌       | ✅      |

:::{dropdown} _A high-level overview of lakehouse technologies._

¹ LaminDB provides snapshot isolation and time travel by managing dataset state as transactional metadata records in Postgres/SQLite rather than mutating existing files. While it does not perform in-place row-level mutations like a SQL database, operations like `Collection.append()` atomically create new collection versions pointing to new, immutable artifacts. This extends core lakehouse ACID guarantees to multimodal datasets without conflict. For more, see {doc}`acid`.

² See the [Developer experience](#time-travel) section for examples.

³ Adding a nullable/optional column without rewriting existing files.

⁴ In LaminDB, via branches (stage, review, merge).

⁵ Raw files have no commit protocol; concurrent writers risk partial writes / last-writer-wins.

⁶ No need for cleaning orphaned files like in Iceberg.

:::

</div>

An increasingly popular approach to addressing Iceberg's limitations is **DuckLake**,[^ducklake-format][^ducklake-v1] developed by the DuckDB team. Rather than storing metadata in files, DuckLake keeps all metadata in a relational database, leaving only parquet files in storage. This gives it cheap writes that can be more frequent, transactions with true concurrent writer support, automatic maintenance via the database's native mechanisms, and native multi-table transactions — all things that are difficult or impossible with Iceberg's file-based metadata. A complementary development in operational workloads is **Lakebase**, which decouples Postgres database compute and storage via Write-Ahead Logs in object storage. This brings serverless, transactional Postgres to live applications and agents, while continuously syncing operational row changes into analytical lakehouses like Delta Lake.

**LaminDB** shares DuckLake's core architectural pattern — using a relational database for metadata and object storage for data — but extends it beyond tables to support any format (parquet, zarr, AnnData, HDF5, images).
This enables unified schema management, data lineage, and registry annotations across complex, multimodal datasets (Table 1).

While Iceberg & DuckLake are based on the parquet format, and LaminDB is format-agnostic, **LanceDB** manages datasets in the Lance format, a columnar format inspired by parquet that's optimized for arrays.[^lancedb-format] To use LanceDB, you need to convert your data into the Lance format.
While LanceDB fits the lakehouse architecture, non-lakehouse architectures for managing array-like data exist, too, in particular, `arraylake` & `tensorstore` for `.zarr` arrays, and `tiledb` for `.tiledb` arrays.[^tiledb] These non-lakehouse technologies are out of scope for this post given the established query engines don't apply to them.

### Decoupled compute & query pushdowns

Rather than locking metadata resolution inside a dedicated query engine or custom SQL driver, LaminDB acts as an independent semantic orchestration layer.

When executing analytical queries, LaminDB first resolves metadata in Postgres or SQLite to yield precise object storage paths. You then pass these paths directly to modern open-source engines like DuckDB, Polars, or PySpark.
Because these engines read native Parquet and Zarr files directly over object storage, query execution retains full optimization benefits:

- Projection Pushdowns: Only downloading requested columns.
- Filter Pushdowns: Reading file footers to execute row-group pruning and skip irrelevant data blocks before fetching them.

This decoupled design ensures that using LaminDB for provenance, lineage, and ACID governance introduces zero performance penalty during data processing and analytics. See {doc}`tables` for implementation details.

### Branching & idempotency

To safely delegate tasks to autonomous agents and distributed teams, data infrastructure must support non-destructive experimentation and repeatable execution.

- **Git-like branching (Write-Audit-Publish):** LaminDB provides database-level branching (`stage`, `review`, `merge`) to support isolated workflows. Agents or developers can perform exploratory writes, schema modifications, or pipeline runs on a dedicated branch without corrupting production data or lineage. Once validated, changes are merged into the main database. For details, see {doc}`manage-changes`.
- **Idempotent execution:** Re-running Python scripts, workflow tasks (Nextflow, Snakemake), or agent traces is inherently safe. LaminDB validates content hashes and metadata prior to writing, ensuring that executing logic multiple times never creates duplicate artifacts or dangling storage objects. For details, see {doc}`idempotency`.

(time-travel)=

### Schema evolution & time travel

To see how these concepts translate into developer experience, let's compare the code required to perform these essential agentic operations—appending data, evolving schemas, and time-traveling.

The first type of write operation we need to perform is adding new data to the system. Rather than just dropping a raw file into a bucket, the following code snippets ensure that a new dataset complies with the schema of the existing dataset, and that it's added in an ACID fashion.

::::::{tab-set}
:::::{tab-item} Iceberg
Atomic and snapshot-isolated. New Parquet files and a snapshot manifest are written to S3; concurrent readers see a consistent state throughout.

```python
table.append(batch)  # batch is a pyarrow dataset
```

:::::

:::::{tab-item} LanceDB
`add()` writes new rows to S3 and automatically increments the table version.

```python
table.add(batch)  # batch is a pyarrow dataset
```

:::::

:::::{tab-item} LaminDB
Atomic and snapshot-isolated. A new parquet file creates a new collection version.

```python
collection.append(batch)  # batch is an artifact
```

:::::
::::::

Similarly, when an analysis requires new features, the following snippets ensure that columns are updated consistently across the entire dataset, and future incoming datasets.

::::::{tab-set}

:::::{tab-item} Iceberg
A new metadata file records the updated schema. Existing Parquet files are not modified; reads of old files return `null` for the new column.

```python
from pyiceberg.types import BooleanType
with table.update_schema() as update:
    update.add_column("QC_PASS", BooleanType())
```

:::::

:::::{tab-item} LanceDB
`add_columns` takes a per-column SQL value expression — hence the `CAST(NULL AS BOOLEAN)` string, which supplies both the value and its type for existing rows.

```python
table.add_columns({"QC_PASS": "CAST(NULL AS BOOLEAN)"})
```

:::::

:::::{tab-item} LaminDB
LaminDB registers the feature in its schema registry, validating all future artifacts instance-wide.

```python
feature = ln.Feature(name="QC_PASS", dtype=bool).save()
collection.schema.add(feature)
```

:::::
::::::

Finally, because agents inevitably make mistakes, we look at how to retrieve a previous version of a dataset via "time travel".

::::::{tab-set}

:::::{tab-item} Iceberg

```python
first_snapshot = table.history()[0].snapshot_id  # access version 0
table.scan(snapshot_id=first_snapshot)
```

:::::

:::::{tab-item} LanceDB

```python
table.checkout(1)             # checkout a previous version
```

:::::

:::::{tab-item} LaminDB

```python
collection.versions.get(version="1")  # get a previous version
```

:::::
::::::

## Database schema & API

LaminDB provides a SQL schema for common metadata entities: {class}`~lamindb.Artifact`, {class}`~lamindb.Collection`, {class}`~lamindb.Transform`, {class}`~lamindb.Feature`, {class}`~lamindb.Record` etc. - see the [API reference](/api) or the [source code](https://github.com/laminlabs/lamindb/tree/main/lamindb/models).

The core metadata schema is extendable through modules, e.g., with basic biological ({class}`~bionty.Gene`, {class}`~bionty.Protein`, {class}`~bionty.CellLine`, etc.) & operational entities (`Biosample`, `Techsample`, `Treatment`, etc.).

Data models are defined in Python using the Django ORM. Django translates them to SQL tables.
[Django](https://github.com/django/django) is one of the most-used & highly-starred projects on GitHub (~1M dependents, ~73k stars) and has been robustly maintained for 15 years.
While the SQLAlchemy ORM has some advantages, Django is the most popular choice for building metadata management systems in the life sciences.

On top of the metadata schema, LaminDB is a Python API that models datasets as artifacts, abstracts storage & database access, data transformations, and ontologies.

## Modules

LaminDB can be extended with modules building on the [Django](https://github.com/django/django) ecosystem. Examples are:

- [bionty](./bionty): Basic biological ontologies, with easy import from >20 public ontologies
- [pertdb](https://github.com/laminlabs/pertdb): Registries for perturbations (compounds, biologics, genetic interventions, etc.)

If you'd like to create your own module:

1. Create a git repository with registries similar to [pertdb](https://github.com/laminlabs/pertdb)
2. Create & deploy migrations via `lamin migrate create` and `lamin migrate deploy`

For more information, see {doc}`docs:setup`.

## Repositories

LaminDB and its plugins consist in open-source Python libraries & publicly hosted metadata assets:

- [lamindb](https://github.com/laminlabs/lamindb): Core library.
- [bionty](https://github.com/laminlabs/bionty): Basic biological ontologies, with easy import from >20 public ontologies
- [pertdb](https://github.com/laminlabs/pertdb): Registries for perturbations (compounds, biologics, genetic interventions, etc.)

Tightly integrated dependencies are available as git submodules [here](https://github.com/laminlabs/lamindb/tree/main/sub), for instance,

- [lamindb-setup](https://github.com/laminlabs/lamindb-setup): Setup & configure LaminDB.
- [lamin-cli](https://github.com/laminlabs/lamin-cli): The CLI.

Use cases / domain-specific repos:

- [lamin-usecases](https://github.com/laminlabs/lamin-usecases): Use cases as visible on the docs.
- [redun-lamin](https://github.com/laminlabs/redun-lamin): Track redun workflow runs with LaminDB.
- [lamin-mlops](https://github.com/laminlabs/lamin-mlops): MLOps use cases (MNIST, W&B, MLflow, Croissant).
- [cellxgene-lamin](https://github.com/laminlabs/cellxgene-lamin): CELLxGENE data and curation.
- [lamin-spatial](https://github.com/laminlabs/lamin-spatial): Spatial data (RxRx, Vitessce).
- [snakemake-lamin](https://github.com/laminlabs/snakemake-lamin): Track Snakemake runs with LaminDB.
- [nf-lamin](https://github.com/laminlabs/nf-lamin): Nextflow integration with LaminDB.

For a comprehensive list of open-sourced software, browse our [GitHub account](https://github.com/laminlabs), for instance,

- [readfcs](https://github.com/laminlabs/readfcs): FCS artifact reader.

There is a public repository for LaminHub:

- [laminhub-public](https://github.com/laminlabs/laminhub-public): Make issues and follow releases of LaminHub, no source code.

## References

[^apache-iceberg]: Apache Software Foundation. Apache Iceberg: The open table format for analytic datasets. [Apache Iceberg](https://iceberg.apache.org/).

[^ducklake-format]: Raasveldt M & Mühleisen H (2025). DuckLake: SQL as a Lakehouse Format. [DuckLake Blog](https://ducklake.select/2025/05/27/ducklake-01/).

[^ducklake-v1]: Raasveldt M & Holanda P (2026). DuckLake v1.0: The Lakehouse Format Built on SQL Reaches Production-Readiness. [DuckLake Blog](https://ducklake.select/2026/04/13/ducklake-10/).

[^delta]: Linux Foundation. Delta Lake: An open-source storage framework that enables building a Lakehouse architecture. [Delta Lake](https://delta.io/).

[^hudi]: Apache Software Foundation. Apache Hudi: Streaming data on data lakes. [Apache Hudi](https://hudi.apache.org/).

[^nessie]: Project Nessie. Nessie: Transactional Catalog for Data Lakes. [Project Nessie](https://projectnessie.org/).

[^databricks]: Databricks (2020). Accurately Building Genomic Cohorts at Scale with Delta Lake and Spark. [Databricks Blog](https://www.databricks.com/blog/2020/09/22/accurately-building-genomic-cohorts-at-scale-with-delta-lake-and-spark.html).

[^lancedb-format]: LanceDB (2024). Lance Format v2.2 Benchmarks: Half the storage, none of the slowdown. [LanceDB Blog](https://lancedb.com/blog/lance-format-v2-2-benchmarks-half-the-storage-none-of-the-slowdown).

[^tiledb]: TileDB (2020). Population Genomics Data with TileDB. [TileDB Blog](https://tiledb.com/blog/population-genomics-data-with-tiledb).
