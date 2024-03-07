**Access data & metadata across storage (files, arrays) & database (SQL) backends.**

- Query & search: {class}`~lamindb.core.Registry.filter`, {class}`~lamindb.core.Registry.search`
- Stage, load & stream artifacts: {class}`~lamindb.Artifact.stage`, {class}`~lamindb.Artifact.load`, {class}`~lamindb.Artifact.backed`
- Manage {class}`~lamindb.Feature`, {class}`~lamindb.FeatureSet`, {class}`~lamindb.ULabel`
- Plug-in custom [schemas](/schemas) & manage schema migrations
- Use array formats in memory & storage: [DataFrame](/tutorial), [AnnData](/data), [MuData](docs:multimodal), [SOMA](docs:cellxgene), ... backed by [parquet](/tutorial), [zarr](/data), [TileDB](docs:cellxgene), [HDF5](/data), [h5ad](/data), [DuckDB](docs:rxrx), ...
- Bridge artifacts and warehousing: {class}`~lamindb.Artifact`, {class}`~lamindb.Collection`
- Leverage out-of-the-box PyTorch data loaders: {meth}`~lamindb.Collection.mapped`
- Version artifacts, collections & transforms

**Track data lineage across notebooks, pipelines & UI: {meth}`~lamindb.track`, {class}`~lamindb.Transform` & {class}`~lamindb.Run`.**

- Execution reports, source code and Python environments for [notebooks & scripts](/track)
- Integrate with workflow managers: [redun](docs:redun), [nextflow](docs:nextflow), [snakemake](docs:snakemake)

**Manage registries for experimental metadata & in-house ontologies, import public ontologies.**

- Use >20 public ontologies with plug-in {mod}`bionty`
- {class}`~bionty.Gene`, {class}`~bionty.Protein`, {class}`~bionty.CellMarker`
- {class}`~bionty.ExperimentalFactor`, {class}`~bionty.CellType`, {class}`~bionty.CellLine`, {class}`~bionty.Tissue`, ...
- Safeguards against typos & duplications

**Validate, standardize & annotate based on registries: {class}`~lamindb.core.CanValidate.validate` & {class}`~lamindb.core.CanValidate.standardize`.**

- Inspect validation failures: {class}`~lamindb.core.CanValidate.inspect`
- Annotate with untyped or typed labels: {class}`~lamindb.core.LabelManager.add`
- Save data & metadata ACID: {class}`~lamindb.Artifact.save`

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
