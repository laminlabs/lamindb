[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![Docs](https://img.shields.io/badge/docs-humans-yellow)](https://docs.lamin.ai)
[![DocsLLMs](https://img.shields.io/badge/docs-LLMs-yellow)](https://docs.lamin.ai/summary.md)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)
[![PyPI Downloads](https://img.shields.io/pepy/dt/lamindb?logo=pypi)](https://pepy.tech/project/lamindb)

# LaminDB - A data framework for biology

Makes your data queryable, traceable, reproducible, and FAIR. One API: lakehouse, lineage, feature store, ontologies, LIMS, ELN.

<details>
<summary>Why?</summary>

Reproducing analytical results or understanding how a dataset or model was created can be a pain.
Training models on historical data, LIMS & ELN systems, orthogonal assays, or datasets from other teams is even harder.
Even maintaining an overview of a project's datasets & analyses is more difficult than it should be.

Biological datasets are typically managed with versioned storage systems, GUI-focused platforms, structureless data lakes, rigid data warehouses (SQL, monolithic arrays), or tabular lakehouses.

LaminDB extends the lakehouse architecture to biological registries & datasets beyond tables (`DataFrame`, `AnnData`, `.zarr`, `.tiledbsoma`, …) with enough structure to enable queries and enough freedom to keep the pace of R&D high.
Moreover, it provides context through data lineage -- tracing data and code, scientists and models -- and abstractions for biological domain knowledge and experimental metadata.

</details>

<img width="800px" src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/BunYmHkyFLITlM5M0005.png">

Highlights:

- **lineage** → track inputs & outputs of notebooks, scripts, functions & pipelines with a single line of code
- **lakehouse** → manage, monitor & validate schemas; query across many datasets
- **feature store** → manage features & labels; leverage batch loading
- **FAIR datasets** → validate & annotate `DataFrame`, `AnnData`, `SpatialData`, `parquet`, `zarr`, …
- **LIMS & ELN** → manage experimental metadata, ontologies & markdown notes
- **unified access** → storage locations (local, S3, GCP, …), SQL databases (Postgres, SQLite) & ontologies
- **reproducible** → auto-version & timestamp execution reports, source code & environments
- **zero lock-in & scalable** → runs in your infrastructure; not a client for a rate-limited REST API
- **integrations** → [vitessce](https://docs.lamin.ai/vitessce), [nextflow](https://docs.lamin.ai/nextflow), [redun](https://docs.lamin.ai/redun), and [more](https://docs.lamin.ai/integrations)
- **extendable** → create custom plug-ins based on the Django ORM

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

Copy [summary.md](https://docs.lamin.ai/summary.md) into an LLM chat and let AI explain or read the [docs](https://docs.lamin.ai).

## Setup

<!-- copied from quick-setup-lamindb.md -->

Install the `lamindb` Python package:

```shell
pip install lamindb
```

Create a LaminDB instance:

```shell
lamin init --modules bionty --storage ./quickstart-data  # or s3://my-bucket, gs://my-bucket
```

Or if you have write access to an instance, connect to it:

```shell
lamin connect account/name
```

## Quickstart

### Lineage

Create a dataset while tracking source code, inputs, outputs, logs, and environment:

```python
import lamindb as ln

ln.track()  # track execution of source code as a run
open("sample.fasta", "w").write(">seq1\nACGT\n")  # create a dataset
ln.Artifact("sample.fasta", key="sample.fasta").save()  # save dataset as an artifact
ln.finish()  # mark the run as finished
```

Running this snippet as a script (`python create-fasta.py`) produces the following data lineage.

```python
artifact = ln.Artifact.get(key="sample.fasta")  # get artifact by key
artifact.view_lineage()
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/EkQATsQL5wqC95Wj0006.png" width="200">

You'll know how that artifact was created and what it's used for. Basic metadata was captured in fields:

```python
artifact.size        # access the size
artifact.created_at  # access the timestamp
# etc.
artifact.describe()  # describe metadata
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/BOTCBgHDAvwglN3U0004.png" width="550">

Here is how to access the content of the artifact:

```python
local_path = artifact.cache()  # return a local path from a cache
object = artifact.load()       # if available for the format, load object into memory
accessor = artifact.open()     # if available for the format, return a streaming accessor
```

And here is how to access its data lineage context:

```python
run = artifact.run                  # get the run record
transform = artifact.run.transform  # get the transform record
```

<details>
<summary>Examples for run & transform.</summary>

```python
run.describe()
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/rJrHr3XaITVS4wVJ0000.png" width="550" />

```python
transform.describe()
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/JYwmHBbgf2MRCfgL0000.png" width="550" />
</details>

### Lake: annotation & queries

You can annotate datasets and samples with features. Let's define some:

```python
from datetime import date

ln.Feature(name="gc_content", dtype=float).save()
ln.Feature(name="experiment_note", dtype=str).save()
ln.Feature(name="experiment_date", dtype=date).save()
```

During annotation, feature names and data types are validated against these definitions:

```python
artifact.features.add_values({
    "gc_content": 0.55,
    "experiment_note": "Looks great",
    "experiment_date": "2025-10-24",
})
```

Now that the data is annotated, you can query for it:

```python
ln.Artifact.filter(experiment_date="2025-10-24").to_dataframe()  # query all artifacts annotated with `experiment_date`
```

You can also query by the metadata that lamindb automatically collects:

```python
ln.Artifact.filter(run=run).to_dataframe()                # query all artifacts created by a run
ln.Artifact.filter(transform=transform).to_dataframe()    # query all artifacts created by a transform
ln.Artifact.filter(size__gt=1e6).to_dataframe()           # query all artifacts bigger than 1MB
```

If you want to include more information into the resulting dataframe, pass `include`.

```python
ln.Artifact.to_dataframe(include="features")  # include the feature annotations
ln.Artifact.to_dataframe(include=["created_by__name", "storage__root"])  # include fields from related registries
```

### Lake ♾️ LIMS ♾️ Sheets

You can create records for the entities underlying your experiments: samples, perturbations, instruments, etc., for example:

```python
sample = ln.Record(name="Sample", is_type=True).save()  # type sample
ln.Record(name="P53mutant1", type=sample).save()        # sample 1
ln.Record(name="P53mutant2", type=sample).save()        # sample 2
```

Define the corresponding features and annotate:

```python
ln.Feature(name="design_sample", dtype=sample).save()
artifact.features.add_values({"design_sample": "P53mutant1"})
```

You can query & search the `Record` registry in the same way as `Artifact` or `Run`.

```python
ln.Record.search("p53").to_dataframe()
```

You can also create relationships of entities and -- if you connect your LaminDB instance to LaminHub -- edit them like Excel sheets in a GUI.

### Lake: versioning

If you change source code or datasets, LaminDB manages their versioning for you.
Assume you run a new version of our `create-fasta.py` script to create a new version of `sample.fasta`.

```python
import lamindb as ln

ln.track()
open("sample.fasta", "w").write(">seq1\nTGCA\n")  # a new sequence
ln.Artifact("sample.fasta", key="sample.fasta", features={"design_sample": "P53mutant1"}).save()  # annotate with the new sample
ln.finish()
```

If you now query by `key`, you'll get the latest version of this artifact.

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

To validate & annotate the content of the dataframe, use a built-in `schema`:

```python
ln.Feature(name="sequence_str", dtype=str).save()  # define a remaining feature
artifact = ln.Artifact.from_dataframe(
    df,
    key="my_datasets/sequences.parquet",
    schema="valid_features"  # validate columns against features
).save()
artifact.describe()
```

Now you know which schema the dataset satisfies. You can filter for datasets by schema and then launch distributed queries and batch loading.

### Lakehouse beyond tables

To validate an `AnnData` with a built-in `schema` call:

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

To validate a `spatialdata` or any other array-like dataset, you need to construct a `Schema`. You can do this by composing the schema of a complicated object from simple `pandera`/`pydantic`-like schemas: [docs.lamin.ai/curate](https://docs.lamin.ai/curate).

### Ontologies

Plugin `bionty` gives you >20 of them as `SQLRecord` registries. This was used to validate the `ENSG` ids in the `adata` just before.

```python
import bionty as bt

bt.CellType.import_source()  # import the default ontology
bt.CellType.to_dataframe()   # your extendable cell type ontology in a simple registry
```

### CLI

Most of the functionality that's available in Python is also available on the command line (and in `R` through `LaminR`). For instance, to upload a file or folder, run:

```shell
lamin save myfile.txt --key examples/myfile.txt
```

### Workflow managers

LaminDB is not a workflow manager, but it integrates well with existing workflow managers and can subsitute them in some settings.

In [github.com/laminlabs/schmidt22](https://github.com/laminlabs/schmidt22) we manage several workflows, scripts, and notebooks to re-construct the project of [Schmidt _el al._ (2022)](https://pubmed.ncbi.nlm.nih.gov/35113687/). A phenotypic CRISPRa screening result is integrated with scRNA-seq data. Here is one of the input artifacts:

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/JvLaK9Icj11eswQn0000.png" width="850">

And here is the lineage of the final result:

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/KQmzmmLOeBN0C8Yk0004.png" width="850">

You can explore it [here](https://lamin.ai/laminlabs/lamindata/artifact/W1AiST5wLrbNEyVq0001).

If you'd like to integrate with Nextflow, Snakemake, or redun, see here: [docs.lamin.ai/pipelines](https://docs.lamin.ai/pipelines)
