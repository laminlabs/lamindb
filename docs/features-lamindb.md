**Access data & metadata across storage (files, arrays) & database (SQL) backends.**

- Query & search: {class}`~lamindb.dev.Registry.filter`, {class}`~lamindb.dev.Registry.search`
- Stage, load or stream files & datasets: {class}`~lamindb.File.stage`, {class}`~lamindb.File.load`, {class}`~lamindb.File.backed`
- Model data using {class}`~lamindb.Feature`, {class}`~lamindb.FeatureSet`, {class}`~lamindb.ULabel`
- Plug-in custom [schemas](/schemas) & manage schema migrations
- Use array formats in memory & storage: [DataFrame](/tutorial), [AnnData](/data), [MuData](docs:multimodal), [SOMA](docs:cellxgene), ... backed by [parquet](/tutorial), [zarr](/data), [TileDB](docs:cellxgene), [HDF5](/data), [h5ad](/data), [DuckDB](docs:rxrx), ...
- Bridge data artifacts and warehousing: {class}`~lamindb.File`, {class}`~lamindb.Dataset`
- Version files, datasets & transforms

**Track data flow across notebooks, pipelines & UI: {meth}`~lamindb.track`, {class}`~lamindb.Transform` & {class}`~lamindb.Run`.**

- Execution reports & source code for [notebooks & scripts](/track)
- Integrate with workflow managers: [redun](docs:redun), [nextflow](docs:nextflow), [snakemake](docs:snakemake)

**Manage registries for experimental metadata & in-house ontologies, import public ontologies.**

- Use >20 public ontologies with plug-in {mod}`lnschema_bionty`
- {class}`~lnschema_bionty.Gene`, {class}`~lnschema_bionty.Protein`, {class}`~lnschema_bionty.CellMarker`
- {class}`~lnschema_bionty.ExperimentalFactor`, {class}`~lnschema_bionty.CellType`, {class}`~lnschema_bionty.CellLine`, {class}`~lnschema_bionty.Tissue`, ...
- Safeguards against typos & duplications

**Validate, standardize & annotate data using registries: {class}`~lamindb.dev.CanValidate.validate` & {class}`~lamindb.dev.CanValidate.standardize`.**

- Inspect validation failures: {class}`~lamindb.dev.CanValidate.inspect`
- Annotate with untyped or typed labels: {class}`~lamindb.dev.LabelManager.add`
- Save data & metadata ACID: {class}`~lamindb.File.save`

**Organize and share data across a mesh of LaminDB instances.**

- Create & load instances like git repos: `lamin init` & `lamin load`
- Zero-copy [transfer](/transfer) data across instances

**Zero lock-in, scalable, auditable, access management, and more.**

- Zero lock-in: LaminDB runs on generic backends server-side and is not a client for "Lamin Cloud"
  - Flexible storage backends (local, S3, GCP, anything [fsspec](https://github.com/fsspec) supports)
  - Currently two SQL backends for managing metadata: SQLite & Postgres
- Scalable: metadata tables support 100s of millions of entries
- Auditable: data & metadata records are hashed, timestamped, and attributed to users (soon to come: LaminDB Log)
- [Access](docs:access) management:
  - High-level access management through Lamin's collaborator roles
  - Fine-grained access management via storage & SQL roles (soon to come: Lamin Vault)
- [Secure](docs:access): embedded in your infrastructure (Lamin has no access to your data & metadata)
- Tested & typed (up to Django Model fields)
- [Idempotent](docs:faq/idempotency) & [ACID](docs:faq/acid) operations
