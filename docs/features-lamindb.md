**Access data & metadata across storage (files, arrays) & database (SQL) backends.**

- Query by & search for anything: {class}`~lamindb.dev.Registry.filter`, {class}`~lamindb.dev.Registry.search`
- Stage, load or stream files & datasets: {class}`~lamindb.File.stage`, {class}`~lamindb.File.load`, {class}`~lamindb.File.backed`
- Model data using built-in concepts: {class}`~lamindb.Feature`, {class}`~lamindb.FeatureSet`, {class}`~lamindb.ULabel`, {class}`~lamindb.Modality`
- Plug-in custom [schemas](/schemas) & manage schema migrations with ease
- Use array formats in memory & storage: [DataFrame](/tutorial), [AnnData](/data), [MuData](docs:multimodal), [SOMA](docs:cellxgene-census), ... backed by [parquet](/tutorial), [zarr](/data), [TileDB](docs:cellxgene-census), [HDF5](/data), [h5ad](/data), [DuckDB](docs:rxrx), ...
- Bridge data artifacts and data warehousing: {class}`~lamindb.File`, {class}`~lamindb.Dataset`

**Track data flow across notebooks, pipelines & UI: {meth}`~lamindb.track`, {class}`~lamindb.Transform` & {class}`~lamindb.Run`.**

- Associate execution reports & source code with [notebooks](/tutorial)
- Integrate with workflow managers: [redun](docs:redun), [nextflow](docs:nextflow), [snakemake](docs:snakemake)

**Manage registries for experimental metadata & in-house ontologies, import public ontologies.**

- Use >20 public ontologies with plug-in {mod}`lnschema_bionty`
- {class}`~lnschema_bionty.Gene`, {class}`~lnschema_bionty.Protein`, {class}`~lnschema_bionty.CellMarker`
- {class}`~lnschema_bionty.ExperimentalFactor`, {class}`~lnschema_bionty.CellType`, {class}`~lnschema_bionty.CellLine`, {class}`~lnschema_bionty.Tissue`, ...
- Safeguards against typos & duplications

**Validate, standardize & annotate data using registries.**

- {class}`~lamindb.dev.CanValidate.validate` & {class}`~lamindb.dev.CanValidate.standardize`, {class}`~lamindb.dev.CanValidate.inspect` validation failures
- annotate with untyped or typed labels: {class}`~lamindb.dev.LabelManager.add`
- save data & metadata ACID: {class}`~lamindb.File.save`

**Organize and share data across a mesh of LaminDB instances.**

- Create & load instances like git repos
- Zero-copy [transfer](/transfer) data across instances

**Zero lock-in, scalable, access management, data versioning, and more.**

- Zero lock-in: LaminDB runs on generic backends server-side and is not a client for "Lamin Cloud"
  - Flexible storage backends (local, S3, GCP, anything [fsspec](https://github.com/fsspec) supports)
  - Currently two SQL backends for managing metadata: SQLite & Postgres
- Scalable: metadata tables support 100s of millions of entries
- Access management:
  - High-level access management through Lamin's collaborator roles
  - Fine-grained access management via embedded storage & SQL roles
- Secure: embedded in your infrastructure (Lamin has no access to your data & metadata)
- File, dataset & transform versioning
- Tested & typed (up to Django Model fields, to come)
- [Idempotent](docs:faq/idempotency) & [ACID](docs:faq/acid) operations
