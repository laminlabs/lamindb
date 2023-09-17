In many organizations, fragmented object stores, SQL databases & ELN/LIMS systems pile up non-findable, inaccessible, hard-to-integrate & non-traceable data.

This also holds for the derived analytical insights, which makes it hard to optimize learning & decision-making across team & ML models.

LaminDB aims to provide a framework that addresses the [key problems](https://lamin.ai/blog/2022/problems) underlying this tendency.

**Key features**

- Unify access to data & metadata across storage & SQL database backends:
  - Query & search by anything: {class}`~lamindb.dev.Registry.filter`, {class}`~lamindb.dev.Registry.search`
  - Stage, load or stream files & datasets: {class}`~lamindb.File.stage`, {class}`~lamindb.File.load`, {class}`~lamindb.File.backed`
  - Model data schema-less or schema-full, add [custom schema plug-ins & manage schema migrations](/schemas)
  - Organize data around learning: {class}`~lamindb.Feature`, {class}`~lamindb.FeatureSet`, {class}`~lamindb.ULabel`, {class}`~lamindb.Modality`
  - Leverage support for common array formats in memory & storage: `DataFrame`, `AnnData`, `MuData`, `pyarrow.Table` backed by `h5ad`, `parquet`, `zarr`, `TileDB`, `HDF5`
  - Bridge immutable data artifacts ({class}`~lamindb.File`) and data warehousing ({class}`~lamindb.Dataset`)
- Track data flow across notebooks, pipelines & UI: {meth}`~lamindb.track`, {class}`~lamindb.Transform` & {class}`~lamindb.Run`
- Manage registries for experimental metadata & ontologies in a simple database:
  - Use >20 public ontologies with plug-in {mod}`lnschema_bionty`
  - For instance, {class}`~lnschema_bionty.Gene`, {class}`~lnschema_bionty.Protein`, {class}`~lnschema_bionty.CellMarker`, {class}`~lnschema_bionty.ExperimentalFactor`, {class}`~lnschema_bionty.CellType` ...
- Validate, standardize & annotate data:
  - {class}`~lamindb.dev.CanValidate.validate` & {class}`~lamindb.dev.CanValidate.standardize`
  - {class}`~lamindb.dev.CanValidate.inspect` validation failures
  - annotate with untyped or typed labels: {class}`~lamindb.dev.LabelManager.add`
  - save data & metadata ACID: {class}`~lamindb.File.save`
- Create DB instances within seconds and share data across a mesh of instances: {mod}`~lamindb.setup`

**For platform builders**

- Zero lock-in: universal backends, LaminDB is not a client for "Lamin Cloud" but runs server-side
  - Flexible storage backends (local, S3, GCP, anything [fsspec](https://github.com/fsspec) supports)
  - Currently two SQL backends for managing metadata: SQLite & Postgres
- Scalable: metadata tables support 100s of millions of entries
- Access management:
  - High-level access management through Lamin's collaborator roles
  - Fine-grained access management via embedded storage & SQL roles
- Secure: embedded in your infrastructure (Lamin has no access to your data & metadata)
- [Idempotent](docs:faq/idempotency) & [ACID](docs:faq/acid) operations
- File, dataset & transform versioning
- Safeguards against typos & duplications when populating registries
