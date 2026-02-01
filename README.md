# LaminDB [![docs](https://img.shields.io/badge/docs-yellow)](https://docs.lamin.ai) [![llms.txt](https://img.shields.io/badge/llms.txt-orange)](https://docs.lamin.ai/llms.txt) [![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb) [![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=PyPI)](https://pypi.org/project/lamindb) [![cran](https://www.r-pkg.org/badges/version/laminr?color=green)](https://cran.r-project.org/package=laminr) [![stars](https://img.shields.io/github/stars/laminlabs/lamindb?style=flat&logo=GitHub&label=&color=gray)](https://github.com/laminlabs/lamindb) [![downloads](https://static.pepy.tech/personalized-badge/lamindb?period=total&units=INTERNATIONAL_SYSTEM&left_color=GRAY&right_color=GRAY&left_text=%E2%AC%87%EF%B8%8F)](https://pepy.tech/project/lamindb)

LaminDB is an open-source data framework for biology to query, trace, and validate datasets and models at scale.
With one API, you get: lakehouse, lineage, feature store, ontologies, bio-registries & formats.

<details>
<summary>Why?</summary>

Reproducing and understanding how datasets, models, and workflows were created is crucial to high-quality R&D — especially as agents increasingly contribute to it.
At the same time, training models across thousands of datasets — from LIMS and ELNs to orthogonal assays and cross-team silos — is now a major learning opportunity, but requires queryable & validated data.
Robustly scaled learning operations need something biology has lacked: an API-first data management framework comparable to git for code or warehouses for tables.

LaminDB fills the gap with a lineage-native data lakehouse that understands bio-registries and formats (`AnnData`, `.zarr`, …).
It provides queries across many datasets with enough freedom to maintain high-paced R&D while automating rich context on top of versioning, change management, and other industry standards.

</details>

<img width="800px" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/BunYmHkyFLITlM5M000C.svg">

Highlights:

- **lineage** → track inputs & outputs of notebooks, scripts, functions & pipelines with a single line of code
- **lakehouse** → manage, monitor & validate schemas; query across many datasets
- **feature store** → manage features & labels; leverage batch loading
- **FAIR datasets** → validate & annotate `DataFrame`, `AnnData`, `SpatialData`, `parquet`, `zarr`, …
- **LIMS & ELN** → manage experimental metadata, ontologies & markdown notes
- **unified access** → single API for storage locations (local, S3, GCP, …), SQL databases (Postgres, SQLite) & ontologies
- **reproducible** → auto-track source code & compute environments with data, code & report versioning
- **zero lock-in** → runs in your infrastructure on open standards (Postgres, SQLite, `parquet`, `zarr`, etc.)
- **scalable** → you hit storage & database directly through your `pydata` or R stack, no REST API involved
- **simple** → just `pip install` a Python package
- **integrations** → [vitessce](https://docs.lamin.ai/vitessce), [nextflow](https://docs.lamin.ai/nextflow), [redun](https://docs.lamin.ai/redun), and [more](https://docs.lamin.ai/integrations)
- **extensible** → create custom plug-ins based on the Django ORM

If you want a GUI: [LaminHub](https://lamin.ai) is a data collaboration hub built on LaminDB similar to how GitHub is built on git.

<details>
<summary>Who uses it?</summary>

Scientists & engineers in pharma, biotech, and academia, including:

- Pfizer – A global BigPharma company with headquarters in the US
- Ensocell Therapeutics – A BioTech with offices in Cambridge, UK, and California
- DZNE – The National Research Center for Neuro-Degenerative Diseases in Germany
- Helmholtz Munich – The National Research Center for Environmental Health in Germany
- scverse – An international non-profit for open-source omics data tools
- The Global Immunological Swarm Learning Network – Research hospitals at U Bonn, Harvard, MIT, Stanford, ETH Zürich, Charite, Mount Sinai, and others

</details>

## Docs

Copy [llms.txt](https://docs.lamin.ai/llms.txt) into an LLM chat and let AI explain or read the [docs](https://docs.lamin.ai).

## Quickstart

Install the Python package:

```shell
pip install lamindb
```

### Query databases

You can browse public databases at [lamin.ai/explore](https://lamin.ai/explore). To query [laminlabs/cellxgene](https://lamin.ai/laminlabs/cellxgene), run:

```python
import lamindb as ln

db = ln.DB("laminlabs/cellxgene")  # a database object for queries
df = db.Artifact.to_dataframe()    # a dataframe listing datasets & models
```

To get a [specific dataset](https://lamin.ai/laminlabs/cellxgene/artifact/BnMwC3KZz0BuKftR), run:

```python
artifact = db.Artifact.get("BnMwC3KZz0BuKftR")  # a metadata object for a dataset
artifact.describe()                             # describe the context of the dataset
```

<details>
<summary>See the output.</summary>
<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/mxlUQiRLMU4Zos6k0001.png" width="550">
</details>

Access the content of the dataset via:

```python
local_path = artifact.cache()  # return a local path from a cache
adata = artifact.load()        # load object into memory
accessor = artifact.open()     # return a streaming accessor
```

You can query by biological entities like `Disease` through plug-in `bionty`:

```python
alzheimers = db.bionty.Disease.get(name="Alzheimer's disease")
df = db.Artifact.filter(diseases=alzheimers).to_dataframe()
```

### Configure your database

You can create a LaminDB instance at [lamin.ai](https://lamin.ai) and invite collaborators.
To connect to a remote instance, run:

```shell
lamin login
lamin connect account/name
```

If you prefer to work with a local SQLite database (no login required), run this instead:

```shell
lamin init --storage ./quickstart-data --modules bionty
```

On the terminal and in a Python session, LaminDB will now auto-connect.

### CLI

To save a file or folder from the command line, run:

```shell
lamin save myfile.txt --key examples/myfile.txt
```

To sync a file into a local cache (artifacts) or development directory (transforms), run:

```shell
lamin load --key examples/myfile.txt
```

Read more: [docs.lamin.ai/cli](https://docs.lamin.ai/cli).

### Lineage: scripts & notebooks

To create a dataset while tracking source code, inputs, outputs, logs, and environment:

```python
import lamindb as ln
# → connected lamindb: account/instance

ln.track()                                              # track code execution
open("sample.fasta", "w").write(">seq1\nACGT\n")        # create dataset
ln.Artifact("sample.fasta", key="sample.fasta").save()  # save dataset
ln.finish()                                             # mark run as finished
```

Running this snippet as a script (`python create-fasta.py`) produces the following data lineage:

```python
artifact = ln.Artifact.get(key="sample.fasta")  # get artifact by key
artifact.describe()      # context of the artifact
artifact.view_lineage()  # fine-grained lineage
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/BOTCBgHDAvwglN3U0004.png" width="550"> <img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/EkQATsQL5wqC95Wj0006.png" width="140">

<details>
<summary>Access run & transform.</summary>

```python
run = artifact.run              # get the run object
transform = artifact.transform  # get the transform object
run.describe()                  # context of the run
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/rJrHr3XaITVS4wVJ0000.png" width="550" />

```python
transform.describe()  # context of the transform
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/JYwmHBbgf2MRCfgL0000.png" width="550" />

</details>

### Lineage: functions & workflows

You can achieve the same traceability for functions & workflows:

```python
import lamindb as ln

@ln.flow()
def create_fasta(fasta_file: str = "sample.fasta"):
    open(fasta_file, "w").write(">seq1\nACGT\n")    # create dataset
    ln.Artifact(fasta_file, key=fasta_file).save()  # save dataset

if __name__ == "__main__":
    create_fasta()
```

Beyond what you get for scripts & notebooks, this automatically tracks function & CLI params and integrates well with established Python workflow managers: [docs.lamin.ai/track](https://docs.lamin.ai/track). To integrate advanced bioinformatics pipeline managers like Nextflow, see [docs.lamin.ai/pipelines](https://docs.lamin.ai/pipelines).

<details>
<summary>A richer example.</summary>

Here is a an automatically generated re-construction of the project of [Schmidt _el al._ (Science, 2022)](https://pubmed.ncbi.nlm.nih.gov/35113687/):

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/KQmzmmLOeBN0C8Yk0004.png" width="850">

A phenotypic CRISPRa screening result is integrated with scRNA-seq data. Here is the result of the screen input:

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/JvLaK9Icj11eswQn0000.png" width="850">

You can explore it [here](https://lamin.ai/laminlabs/lamindata/artifact/W1AiST5wLrbNEyVq) on LaminHub or [here](https://github.com/laminlabs/schmidt22) on GitHub.

</details>

### Labeling & queries by fields

You can label an artifact by running:

```python
my_label = ln.ULabel(name="My label").save()   # a universal label
project = ln.Project(name="My project").save() # a project label
artifact.ulabels.add(my_label)
artifact.projects.add(project)
```

Query for it:

```python
ln.Artifact.filter(ulabels=my_label, projects=project).to_dataframe()
```

You can also query by the metadata that lamindb automatically collects:

```python
ln.Artifact.filter(run=run).to_dataframe()              # by creating run
ln.Artifact.filter(transform=transform).to_dataframe()  # by creating transform
ln.Artifact.filter(size__gt=1e6).to_dataframe()         # size greater than 1MB
```

If you want to include more information into the resulting dataframe, pass `include`.

```python
ln.Artifact.to_dataframe(include=["created_by__name", "storage__root"])  # include fields from related registries
```

Note: The query syntax for `DB` objects and for your default database is the same.

### Queries by features

You can annotate datasets and samples with features. Let's define some:

```python
from datetime import date

ln.Feature(name="gc_content", dtype=float).save()
ln.Feature(name="experiment_note", dtype=str).save()
ln.Feature(name="experiment_date", dtype=date, coerce=True).save()  # accept date strings
```

During annotation, feature names and data types are validated against these definitions:

```python
artifact.features.add_values({
    "gc_content": 0.55,
    "experiment_note": "Looks great",
    "experiment_date": "2025-10-24",
})
```

Query for it:

```python
ln.Artifact.filter(experiment_date="2025-10-24").to_dataframe()  # query all artifacts annotated with `experiment_date`
```

If you want to include the feature values into the dataframe, pass `include`.

```python
ln.Artifact.to_dataframe(include="features")  # include the feature annotations
```

### Lake ♾️ LIMS ♾️ Sheets

You can create records for the entities underlying your experiments: samples, perturbations, instruments, etc., for example:

```python
sample = ln.Record(name="Sample", is_type=True).save()  # create entity type: Sample
ln.Record(name="P53mutant1", type=sample).save()        # sample 1
ln.Record(name="P53mutant2", type=sample).save()        # sample 2
```

Define features and annotate an artifact with a sample:

```python
ln.Feature(name="design_sample", dtype=sample).save()
artifact.features.add_values({"design_sample": "P53mutant1"})
```

You can query & search the `Record` registry in the same way as `Artifact` or `Run`.

```python
ln.Record.search("p53").to_dataframe()
```

You can also create relationships of entities and edit them like Excel sheets in a GUI via LaminHub.

### Data versioning

If you change source code or datasets, LaminDB manages versioning for you.
Assume you run a new version of our `create-fasta.py` script to create a new version of `sample.fasta`.

```python
import lamindb as ln

ln.track()
open("sample.fasta", "w").write(">seq1\nTGCA\n")  # a new sequence
ln.Artifact("sample.fasta", key="sample.fasta", features={"design_sample": "P53mutant1"}).save()  # annotate with the new sample
ln.finish()
```

If you now query by `key`, you'll get the latest version of this artifact with the latest version of the source code linked with previous versions of artifact and source code are easily queryable:

```python
artifact = ln.Artifact.get(key="sample.fasta")  # get artifact by key
artifact.versions.to_dataframe()                # see all versions of that artifact
```

### Lakehouse ♾️ feature store

Here is how you ingest a `DataFrame`:

```python
import pandas as pd

df = pd.DataFrame({
    "sequence_str": ["ACGT", "TGCA"],
    "gc_content": [0.55, 0.54],
    "experiment_note": ["Looks great", "Ok"],
    "experiment_date": [date(2025, 10, 24), date(2025, 10, 25)],
})
ln.Artifact.from_dataframe(df, key="my_datasets/sequences.parquet").save()  # no validation
```

To validate & annotate the content of the dataframe, use the built-in schema `valid_features`:

```python
ln.Feature(name="sequence_str", dtype=str).save()  # define a remaining feature
artifact = ln.Artifact.from_dataframe(
    df,
    key="my_datasets/sequences.parquet",
    schema="valid_features"  # validate columns against features
).save()
artifact.describe()
```

You can filter for datasets by schema and then launch distributed queries and batch loading.

### Lakehouse beyond tables

To validate an `AnnData` with built-in schema `ensembl_gene_ids_and_valid_features_in_obs`, call:

```python
import anndata as ad
import numpy as np

adata = ad.AnnData(
    X=pd.DataFrame([[1]*10]*21).values,
    obs=pd.DataFrame({'cell_type_by_model': ['T cell', 'B cell', 'NK cell'] * 7}),
    var=pd.DataFrame(index=[f'ENSG{i:011d}' for i in range(10)])
)
artifact = ln.Artifact.from_anndata(
    adata,
    key="my_datasets/scrna.h5ad",
    schema="ensembl_gene_ids_and_valid_features_in_obs"
)
artifact.describe()
```

To validate a `spatialdata` or any other array-like dataset, you need to construct a `Schema`. You can do this by composing simple `pandera`-style schemas: [docs.lamin.ai/curate](https://docs.lamin.ai/curate).

### Ontologies

Plugin `bionty` gives you >20 public ontologies as `SQLRecord` registries. This was used to validate the `ENSG` ids in the `adata` just before.

```python
import bionty as bt

bt.CellType.import_source()  # import the default ontology
bt.CellType.to_dataframe()   # your extendable cell type ontology in a simple registry
```

Read more: [docs.lamin.ai/manage-ontologies](https://docs.lamin.ai/manage-ontologies).
