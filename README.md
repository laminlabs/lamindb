[![docs](https://img.shields.io/badge/docs-yellow)](https://docs.lamin.ai) [![llms.txt](https://img.shields.io/badge/llms.txt-orange)](https://docs.lamin.ai/llms.txt) [![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb) [![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=PyPI)](https://pypi.org/project/lamindb) [![cran](https://www.r-pkg.org/badges/version/laminr?color=green)](https://cran.r-project.org/package=laminr) [![stars](https://img.shields.io/github/stars/laminlabs/lamindb?style=flat&logo=GitHub&label=&color=gray)](https://github.com/laminlabs/lamindb) [![downloads](https://static.pepy.tech/personalized-badge/lamindb?period=total&units=INTERNATIONAL_SYSTEM&left_color=GRAY&right_color=GRAY&left_text=%E2%AC%87%EF%B8%8F)](https://pepy.tech/project/lamindb)

# LaminDB: Data management for traceable, multimodal AI

LaminDB is an open-source data management tool that makes it easy to query, trace and govern datasets across diverse storage formats and locations.
Like git, LaminDB runs anywhere and captures relevant context about your work.
This includes the data flow through models and analyses, the entities and notes defining your work, and the features & schemas of datasets.
It takes about a minute to install LaminDB and create a database on your laptop or in the cloud.

<details>
<summary>Why?</summary>

1. Untraceable results cannot be trusted, especially when non-verifiable tasks are delegated to agents.
2. Without effective access to multimodal data, models burn tokens or fail entirely.
3. Without governing changes to data akin to governing changes to software with git, it's hard to evaluate agents, debug their mistakes, and safely merge their contributions.

Especially in life sciences, hard-to-verify tasks are abundant, data formats are heterogeneous, and teams need end-to-end traceability for GxP compliance (21 CFR Part 11 and EU Annex 11).

Existing data lakehouses aren't a comprehensive solution. They only address the 2nd and 3rd of these problems for a single table or array format, such as `parquet` or `zarr`.

</details>

<img width="800px" alt="lamindb-schematic" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/BunYmHkyFLITlM5M000E.svg">

How?

- **lineage** → trace results across agent sessions, notebooks, scripts & workflows
- **lakehouse** → manage datasets in any format (`parquet`, `zarr`, ...) with time travel, schema evolution & [ACID guarantees](https://docs.lamin.ai/acid); load/query with your favorite engine (Polars, DuckDB, `annbatch`, ...)
- **LIMS & ELN** → unified schema-based records management with support for ontologies & notes
- **FAIR datasets** → validate & annotate files, `DataFrame`, `AnnData`, `SpatialData`, …
- **governance** → [manage changes](https://docs.lamin.ai/manage-changes) via branching & versioning

Architecture? → [docs](https://docs.lamin.ai/architecture)

- **zero lock-in** → open source, using open standards (metadata in SQLite/Postgres, data in `parquet`, `zarr`, etc.)
- **scalable** → hit storage & database directly through your `pydata` or R stack, no REST API involved
- **simple** → `pip install lamindb` or `install.packages('laminr')` · no Docker required, no separate backend
- **unified access** → local, S3, GCP, etc. · Postgres, SQLite · ontologies
- **distributed** → [zero-copy data sharing across databases & storage](https://docs.lamin.ai/transfer)
- **reproducible** → [track](https://docs.lamin.ai/track) agent traces, source code & compute environments
- **ACID** → snapshot isolation & time travel via transactional metadata records across datasets in any format (`parquet`, `zarr`, etc.)
- **idempotent** → [re-run](https://docs.lamin.ai/idempotency) logic without worries about duplications or overwrites
- **decoupled compute** → run your favorite engine (Polars, DuckDB, data loaders, ...) with all its benefits
- **integrations** → [bio ontologies](https://docs.lamin.ai/bionty), [git](https://docs.lamin.ai/track#sync-code-with-git), [nextflow](https://docs.lamin.ai/nextflow), [vitessce](https://docs.lamin.ai/vitessce), [redun](https://docs.lamin.ai/redun), and [more](https://docs.lamin.ai/integrations)
- **extensible** → create custom plug-ins based on the Django ORM, the basis for LaminDB's registries

<details>
<summary>Who?</summary>

Scientists and engineers at leading research institutions and biotech companies, including:

- **Industry** → Pfizer, Altos Labs, RefinedScience, Ensocell Therapeutics, ...
- **Academia & Research** → scverse, DZNE (National Research Center for Neuro-Degenerative Diseases), Helmholtz Munich (National Research Center for Environmental Health), ...
- **Research Hospitals** → Global Immunological Swarm Learning Network: Harvard, MIT, Stanford, ETH Zürich, Charité, U Bonn, Mount Sinai, ...

From personal research projects to pharma-scale deployments managing petabytes of data across:

entities | OOMs
--- | ---
observations & datasets | 10¹² & 10⁶
runs & transforms| 10⁹ & 10⁵
proteins & genes | 10⁹ & 10⁶
biosamples & species | 10⁵ & 10²
... | ...

</details>

UI, permissions, audit logs? [LaminHub](https://lamin.ai) is a collaboration hub built on LaminDB similar to how GitHub is built on git.

## Quickstart

To install the Python package with recommended dependencies, use:

```shell
pip install lamindb
```

<details>
<summary>Install with minimal dependencies.</summary>

The `lamindb` package adds data-science related dependencies through the `[full]` extra, see [here](https://github.com/laminlabs/lamindb/blob/2cc91adcf6077c5af69c1a098699085bb0844083/pyproject.toml#L30-L49).

For a minimal install of the `lamindb` namespace, use:

```shell
pip install lamindb-core
```

</details>

Agent? See [`.agents/`](https://github.com/laminlabs/lamindb/tree/main/lamindb/.agents), [`docs/`](https://github.com/laminlabs/lamindb/tree/main/docs) and [llms.txt](https://docs.lamin.ai/llms.txt).

### Query databases

You can explore public databases [here](https://lamin.ai/explore). To query [laminlabs/cellxgene](https://lamin.ai/laminlabs/cellxgene), run:

```python
import lamindb as ln

db = ln.DB("laminlabs/cellxgene")  # a database object for queries
```

Get a **single artifact** in `h5ad` format:

```python
artifact = db.Artifact.get("x4xsAnj6Uwbq7Bhq")  # a metadata object for an artifact
artifact.describe()                             # describe the metadata
```

<details>
<summary>See the output.</summary>

Features are validated and annotated based on ontologies.

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/mxlUQiRLMU4Zos6k0002.png" width="600">

To explore the dataset on the UI, go [here](https://lamin.ai/laminlabs/cellxgene/artifact/x4xsAnj6Uwbq7Bhq).
</details>

Access the content of the artifact via one of:

```python
local_path = artifact.cache()  # return a local path from a cache
adata = artifact.load()        # load object into memory
accessor = artifact.open()     # return a streaming accessor
```

Query the set of artifacts that's annotated by `liver`:

```python
db.Artifact.filter(tissues__name="liver").to_dataframe()
```

Get a **versioned collection** of artifacts in `zarr` format:

```python
db = ln.DB("laminlabs/arrayloader-benchmarks")
collection = db.Collection.get(key="Tahoe100M_zarr")
collection.artifacts.to_dataframe()
```

<details>
<summary>Train a model on a sharded dataset.</summary>

This example uses `annbatch` with the pre-shuffled sharded `Tahoe100M_zarr` collection ([explore](https://lamin.ai/laminlabs/arrayloader-benchmarks/collection/LaJOdLd0xZ3v5ZBw)):

<!-- #region -->

```python
import anndata as ad
import zarr
from annbatch import Loader

paths = [artifact.cache() for artifact in collection.artifacts.all()]
loader = Loader(shuffle=True, batch_size=4096, chunk_size=256, preload_nchunks=64)
loader.add_datasets(
    datasets=[ad.io.sparse_dataset(zarr.open(p)["X"]) for p in paths],
    obs=[ad.io.read_elem(zarr.open(p)["obs"]) for p in paths],
)
for batch in loader:
    pass
```

<!-- #endregion -->

</details>

Get a **validated, versioned collection** of artifacts in `parquet` format:

```python
db = ln.DB("laminlabs/1000genomes")
collection = db.Collection.get("hVu9puwdRGskm1I6")
collection.schema.describe()
```

<details>
<summary>Ad-hoc query a sharded dataset.</summary>

See the validating schema of the collection ([explore](https://lamin.ai/laminlabs/1000genomes/collection/hVu9puwdRGskm1I6)):

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/5xy5hLxrhA0Bob3E0001.png" width="550">

This example uses `polars` to filter the collection of parquet files:

<!-- #region -->

```python
chrom = "1"
lo, hi = 150_000_000, 200_000_000

with collection.open(engine="polars") as df:
    filtered = df.filter(
        (pl.col("chrom") == chrom) & (pl.col("pos") >= lo) & (pl.col("pos") <= hi)
    ).collect()
```

<!-- #endregion -->

</details>

### Create artifacts

To initialize a new database in a development directory, run:

```shell
lamin init --modules bionty  # mount module bionty
```

To connect to an existing database - you can create a hosted one at [lamin.ai](https://lamin.ai) - run:

```shell
lamin login
lamin connect account/name  # tip: add flag `--here` to scope to current directory
```

To save a file or folder:

```python
import lamindb as ln
# → connected lamindb: account/instance

open("sample.fasta", "w").write(">seq1\nACGT\n")        # create dataset
ln.Artifact("sample.fasta", key="sample.fasta").save()  # save dataset
```

To save a file or folder via the CLI, run:

```shell
lamin save sample.fasta --key sample.fasta
```

:::{dropdown} Load it.

To load an artifact via the CLI into a local cache, run:

```shell
lamin load --key sample.fasta
```

Read more about the CLI: [docs.lamin.ai/cli](https://docs.lamin.ai/cli).

:::

### Trace data, code & agents

To track **agent sessions**, install the [library skill](https://github.com/laminlabs/lamindb/tree/main/lamindb/.agents): `uvx library-skills --skill lamindb`.

To track source code, inputs, outputs, logs, and environment of a **script or notebook**, use `ln.track()` and `ln.finish()`:

```python
import lamindb as ln
# → connected lamindb: account/instance

ln.track()                                              # track code execution
open("sample.fasta", "w").write(">seq1\nACGT\n")        # create dataset
ln.Artifact("sample.fasta", key="sample.fasta").save()  # save dataset
ln.finish()                                             # mark run as finished
```

If running this snippet as a script (`python create_fasta.py`), you'll get the following data lineage:

```python
artifact = ln.Artifact.get(key="sample.fasta")  # get artifact by key
artifact.describe()      # context of the artifact
artifact.view_lineage()  # fine-grained lineage
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/BOTCBgHDAvwglN3U0004.png" width="550"> <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/EkQATsQL5wqC95Wj0006.png" width="140">

Watch a mini video: [youtu.be/yK3ODFZLL1A](https://youtu.be/yK3ODFZLL1A)

<details>
<summary>Access run & transform.</summary>

```python
run = artifact.run                     # get the run object
transform = artifact.run.transform     # get the transform object
run.describe()                  # context of the run
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/rJrHr3XaITVS4wVJ0000.png" width="550" />

```python
transform.describe()  # context of the transform
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/JYwmHBbgf2MRCfgL0000.png" width="550" />

</details>

To track a **workflow**, decorate a function:

<!-- #skip_laminr -->

```python
import lamindb as ln

@ln.flow()
def create_fasta(fasta_file: str = "sample.fasta"):
    open(fasta_file, "w").write(">seq1\nACGT\n")    # create dataset
    ln.Artifact(fasta_file, key=fasta_file).save()  # save dataset

if __name__ == "__main__":
    create_fasta()
```

<!-- #end_skip_laminr -->

Beyond what you get for scripts & notebooks, this automatically tracks function & CLI params and integrates well with established Python workflow managers: [docs.lamin.ai/track](https://docs.lamin.ai/track). To integrate advanced bioinformatics pipeline managers like Nextflow, see [docs.lamin.ai/pipelines](https://docs.lamin.ai/pipelines).

<details>
<summary>A richer example.</summary>

Here is a re-construction of the project of [Schmidt _et al._ (Science, 2022)](https://pubmed.ncbi.nlm.nih.gov/35113687/):

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/KQmzmmLOeBN0C8Yk0004.png" width="850">

A phenotypic CRISPRa screening result is integrated with scRNA-seq data. Here is the result of the screen input:

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/JvLaK9Icj11eswQn0000.png" width="850">

You can explore it [here](https://lamin.ai/laminlabs/lamindata/artifact/W1AiST5wLrbNEyVq).

</details>

### Annotate & query artifacts

To annotate and query by **entities**:

```python
# create entities
my_label = ln.ULabel(name="My universal label").save()
my_project = ln.Project(name="My project label").save()

# annotate
artifact.ulabels.add(my_label)
artifact.projects.add(my_project)

# query
ln.Artifact.filter(ulabels=my_label, projects=my_project).to_dataframe()
```

This also works for custom entities, for example, ontology terms:

```python
import bionty as bt

# import an entity from a public ontology source, e.g., by name or an ontology ID
hek293 = bt.CellLine.from_source(name="HEK293").save()
hras = bt.Gene.from_source(ensembl_gene_id="ENSG00000174775").save()

# annotate
artifact.cell_lines.add(hek293)
artifact.genes.add(hras)

# query
ln.Artifact.filter(cell_lines=hek293, genes=hras).to_dataframe(
    include=["cell_lines__name", "genes__symbol"]
)
```

To annotate & query by **features**:

```python
from datetime import date

# define features
gc_content = ln.Feature(name="gc_content", dtype=float).save()
seq_date = ln.Feature(name="seq_date", dtype=date).save()

# annotate
artifact.features.set_values({
    gc_content: 0.55,
    seq_date: date.today(),
})

# query
ln.Artifact.filter(seq_date == date.today()).to_dataframe(include="features")
```

To query by **auto-tracked metadata**:

```python
# by storage location and size
ln.Artifact.filter(storage__root=ln.settings.storage.root)
ln.Artifact.filter(size__gt=1e6).to_dataframe()

# by data lineage
ln.Artifact.filter(run=run).to_dataframe()
ln.Artifact.filter(run__transform=transform).to_dataframe()

# include more information in resulting dataframe
ln.Artifact.to_dataframe(include=["run__transform__key", "storage__root"])
```

Read more: [docs.lamin.ai/organize](https://docs.lamin.ai/organize) or [docs.lamin.ai/query-search](https://docs.lamin.ai/query-search)

### Create records & notes

You can create **record types** to model metadata for any entity:

```python
# create record types
experiment_type = ln.Record(name="Experiments", is_type=True).save()
sample_type = ln.Record(name="Samples", is_type=True).save()

# create an experiment
experiment1 = ln.Record(
    name="Experiment 1",
    features={seq_date: date.today()},
    type=experiment_type
).save()

# create a sample that links to "Experiment 1" via a feature
experiment = ln.Feature(name="experiment", dtype=experiment_type).save()
sample1 = ln.Record(
    features={gc_content: 0.5, experiment: experiment1},
    type=sample_type
).save()

# export all experiments and all samples
experiment_type.to_dataframe()
sample_type.to_dataframe()
```

You can save **markdown notes** as records:

```shell
lamin save <topic>/<my-note.md>
```

To load a record back as a markdown file, run:

```shell
lamin load https://lamin.ai/laminlabs/lamindata/record/mNDJgWFrkWQVW3ox
```

Watch a mini video: [youtu.be/NRzVQXJaRH8](https://youtu.be/NRzVQXJaRH8)

Read more: [docs.lamin.ai/manage-records](https://docs.lamin.ai/manage-records).

### Validate datasets

Here is how you ingest a `DataFrame`:

```python
import pandas as pd

df = pd.DataFrame({
    "seq_str": ["ACGT", "TGCA"],
    "gc_content": [0.55, 0.54],
    "experiment_note": ["Looks great", "Ok"],
    "seq_date": [date(2025, 10, 24), date(2025, 10, 25)],
})
ln.Artifact.from_dataframe(df, key="my_datasets/sequences.parquet").save()  # no validation
```

The simplest way to validate & annotate the content of the dataframe is the built-in schema `valid_features`:

```python
ln.Feature(name="seq_str", dtype=str).save()  # define remaining features
ln.Feature(name="experiment_note", dtype=str).save()
artifact = ln.Artifact.from_dataframe(
    df,
    key="my_datasets/sequences.parquet",
    schema="valid_features"  # validate columns against features
).save()
artifact.describe()
```

You can filter for datasets by schema and then launch distributed queries or batch load distributed datasets. For tables, see: [docs.lamin.ai/tables](https://docs.lamin.ai/tables). For arrays, see: [docs.lamin.ai/arrays](https://docs.lamin.ai/arrays).

Watch a mini video: [youtu.be/Ji6E7hTnReQ](https://youtu.be/Ji6E7hTnReQ)

To validate an `AnnData`, call:

```python
import anndata as ad
import numpy as np
import pandas as pd

adata = ad.AnnData(
    X=np.ones((21, 10)),
    obs=pd.DataFrame({'cell_type_by_model': ['T cell', 'B cell', 'NK cell'] * 7}),
    var=pd.DataFrame(index=[f'ENSG{i:011d}' for i in range(10)])
)
artifact = ln.Artifact.from_anndata(
    adata,
    key="my_datasets/scrna.h5ad",
    schema="ensembl_gene_ids_and_valid_features_in_obs"
).save()
artifact.describe()
```

To validate a `SpatialData` or any other array-like dataset, you need to construct a `Schema`. You can do this by composing simple `pandera`-style schemas: [docs.lamin.ai/curate](https://docs.lamin.ai/curate).

### Use branching & versioning

LaminDB co-versions code and datasets for you.
If you edit and run the `create_fasta.py` script, you'll automatically create a new version of the transform and the `sample.fasta` artifact.

<details>
<summary>The edited script</summary>

```python
# create_fasta.py
import lamindb as ln

ln.track()
open("sample.fasta", "w").write(">seq1\nTGCA\n")  # a new sequence
ln.Artifact("sample.fasta", key="sample.fasta", features={"experiment": "Experiment 1"}).save()  # annotate with the new experiment
ln.finish()
```
</details>

```python
artifact_latest = ln.Artifact.get(key="sample.fasta")  # pass version for a previous version: ln.Artifact.get(key="sample.fasta", version="1.0")
artifact_latest.versions.to_dataframe()                # all versions of that artifact
artifact_latest.run.transform.versions.to_dataframe()  # all versions of the transform that created the artifact
```

To isolate changes, create a contribution branch and switch to it as in `git`:

```shell
lamin switch -c my_branch
```

To merge a contribution branch into `main`, run:

```shell
lamin switch main  # switch to the main branch
lamin merge my_branch  # merge contribution branch into main
```

Read more: [docs.lamin.ai/manage-changes](https://docs.lamin.ai/manage-changes).

Watch a mini video: [youtu.be/rzRwcMj6-fc](https://youtu.be/rzRwcMj6-fc)

### Share datasets

To share data in a lineage-aware way, transfer objects from a source database to your default database:

```python
db = ln.DB("laminlabs/lamindata")
artifact = db.Artifact.get(key="example_datasets/mini_immuno/dataset1.h5ad")
artifact.save()
```

This is zero-copy for the artifact's data in storage. Read more: [docs.lamin.ai/transfer](https://docs.lamin.ai/transfer).

### Manage ontologies

`bionty` gives you >20 public ontologies so that you don't have build ontologies for basic biological entities from scratch. The `Gene` registry was used to validate the `ENSG` ids in the `adata` just before. Here is how you import a public ontology source into the `CellType` registry:

```python
import bionty as bt

bt.CellType.import_source()  # import the default ontology
bt.CellType.to_dataframe()   # your extensible cell type ontology in a simple registry
```

You can then create objects, e.g. for labeling, analogous to `ULabel` or `Project`:

```python
t_cell = bt.CellType.get(name="T cell")
artifact.cell_types.add(t_cell)
```

Read more: [docs.lamin.ai/manage-ontologies](https://docs.lamin.ai/manage-ontologies).

Watch a mini video: [youtu.be/3vpWjHj3Kw8](https://youtu.be/3vpWjHj3Kw8)
