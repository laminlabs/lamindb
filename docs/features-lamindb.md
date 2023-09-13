**Why?**

In many organizations, fragmented object stores, SQL databases & ELN/LIMS systems pile up non-findable, inaccessible & hard-to-integrate data.

This holds for key analytical insights - just another data representation - making it hard to optimize learning & decision-making.

LaminDB attempts to provide a unified framework that address the [key problems](https://lamin.ai/blog/2022/problems) underlying this tendency.

**What?**

- Unified API to manage data in storage & metadata in SQL database
  - Access: query & search, stage & stream
  - Validate & standardize metadata of vectors & arrays
  - Register & annotate: use `.labels.add()` to annotate with untyped or typed labels
- Track data flow (provenance/ lineage) across notebooks, pipelines & UI
- Manage ontologies, embed custom ontologies in public knowledge
- Model data schema-less or schema-full, add custom schema plug-ins, manage schema migrations
- No silos: create DB instances within seconds and share data across a mesh of instances
- Awareness of common array formats in memory & storage: `DataFrame`, `AnnData`, `MuData`, `pyarrow.Table` backed by `parquet`, `zarr`, `TileDB`, `HDF5`
  - `.backed()` provides access to connectors that allow stream-based queries
  - `.join()` allows auto-concatenating DataFrames

**How?**

- Zero lock-in: universal, simple storage formats, LaminDB is **not** a client for "Lamin Cloud" but can run server-side allowing you to build your own apps
- Scalable: metadata tables support 100s of millions of entries at reasonable query speed, otherwise Lamin relies on object stores
- Flexible storage backends (local, S3, GCP, anything `fsspec` supports)
- Currently two SQL backends: SQLite & Postgres (more to come)
- Fine-grained access management via embedded storage & SQL restrictions & high-level access management through Lamin's collaborator roles
- Secure: fully secure & embedded in your infrastructure & Lamin has no access to your data & metadata
- Idempotent & ACID operations
- File & dataset versioning
- Numerous safeguards against typos & duplications when populating
- A small number of high-quality dependencies
