[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![Docs](https://img.shields.io/badge/docs-humans-yellow)](https://docs.lamin.ai)
[![DocsLLMs](https://img.shields.io/badge/docs-LLMs-yellow)](https://docs.lamin.ai/summary.md)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)
[![PyPI Downloads](https://img.shields.io/pepy/dt/lamindb?logo=pypi)](https://pepy.tech/project/lamindb)

# LaminDB - A data framework for biology

Leverage queryable, traceable, reproducible, and [FAIR](https://en.wikipedia.org/wiki/FAIR_data) data. All in one: lakehouse, feature store, ontologies, lineage, LIMS, ELN.

<details>
<summary>Why?</summary>

Reproducing analytical results or understanding how a dataset or model was created can be a pain.
Training models on historical data, LIMS & ELN systems, orthogonal assays, or datasets from other teams is even harder.
Even maintaining an overview of a project's datasets & analyses is more difficult than it should be.

Biological datasets are typically managed with versioned storage systems, GUI-focused platforms, structureless data lakes, rigid data warehouses (SQL, monolithic arrays), or tabular lakehouses.

LaminDB extends the lakehouse architecture to biological registries & datasets beyond tables (`DataFrame`, `AnnData`, `.zarr`, `.tiledbsoma`, ...) with enough structure to enable queries and enough freedom to keep the pace of R&D high.
Morever, it provides context through data lineage -- tracing data and code, scientists and models -- and abstractions for biological domain knowledge and experimental metadata.

</details>

**Highlights.**

- **lineage:** track inputs & outputs of notebooks, scripts, functions & pipelines with a single line of code
- **lakehouse**: manage, monitor & validate schemas; query across many datasets
- **feature store**: manage features & labels; leverage batch loading
- **FAIR bio datasets**: validate & annotate `DataFrame`, `AnnData`, `SpatialData`, ... backed by `parquet`, `zarr`, HDF5, ...
- **LIMS & ELN**: manage experimental metadata, ontologies & markdown documentation based on Django the ORM
- **unified access:** storage locations (local, S3, GCP, ...), SQL databases (Postgres, SQLite) & ontologies
- **reproducible & auditable:** auto-version & timestamp execution reports, source code & environments; attribute records to users
- **integrations:** [vitessce](https://docs.lamin.ai/vitessce), [nextflow](https://docs.lamin.ai/nextflow), [redun](https://docs.lamin.ai/redun), and [more](https://docs.lamin.ai/integrations)
- **zero lock-in & scalable:** runs in your infrastructure; not a client for a rate-limited REST API
- **extendable:** create custom plug-ins based on the Django ecosystem
- **production-ready:** used in BigPharma, BioTech, hospitals & top labs

If you want a GUI, connect your LaminDB instance to LaminHub and close the drylab-wetlab feedback loop: [lamin.ai](https://lamin.ai).

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
lamin init --storage ./quickstart-data  # or s3://my-bucket, gs://my-bucket
```

Or if you have write access to an instance, connect to it:

```shell
lamin connect account/name
```

## Quickstart

<!-- copied from preface.md -->

Track a script or notebook run with source code, inputs, outputs, logs, and environment.

<!-- copied from py-quickstart.py -->

```python
import lamindb as ln

ln.track()  # track a run
open("sample.fasta", "w").write(">seq1\nACGT\n")
ln.Artifact("sample.fasta", key="sample.fasta").save()  # create an artifact
ln.finish()  # finish the run
```

<!-- from here on, slight deviation from preface.md, where all this is treated in the walk through in more depth -->

This code snippet creates an artifact, which can store a dataset or model as a file or folder in various formats.
Running the snippet as a script (`python create-fasta.py`) produces the following data lineage.

```python
artifact = ln.Artifact.get(key="sample.fasta")  # query artifact by key
artifact.view_lineage()
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/EkQATsQL5wqC95Wj0005.png" width="250">

You'll know how that artifact was created and what it's used for ([interactive visualization](https://lamin.ai/laminlabs/lamindata/artifact/8incOOgjn6F0K1TS)) in addition to capturing basic metadata:

```python
artifact.describe()
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/BOTCBgHDAvwglN3U0002.png" width="550">

You can organize datasets with validation & annotation of any kind of metadata to then access them via queries & search. Here is a more [comprehensive example](https://lamin.ai/laminlabs/lamindata/artifact/9K1dteZ6Qx0EXK8g):

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/6sofuDVvTANB0f480002.png" width="850">

To annotate an artifact with a label, use:

```python
my_experiment = ln.Record(name="My experiment").save()  # create a label in the universal label ontology
artifact.records.add(my_experiment)  # annotate the artifact with the label
```

To query for a set of artifacts, use the `filter()` statement.

```python
ln.Artifact.filter(records=my_experiment, suffix=".fasta").to_dataframe()  # query by suffix and the ulabel we just created
ln.Artifact.filter(transform__key="create-fasta.py").to_dataframe()  # query by the name of the script we just ran
```

If you have a structured dataset like a `DataFrame`, an `AnnData`, or another array, you can validate the content of the dataset (and parse annotations).
Here is [an example for a dataframe](https://docs.lamin.ai/tutorial#validate-an-artifact).

With a large body of validated datasets, you can then access data through distributed queries & batch streaming, see here: [docs.lamin.ai/arrays](https://docs.lamin.ai/arrays).
