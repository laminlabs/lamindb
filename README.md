[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB: Manage R&D data & analyses

_Curate, store, track, query, integrate, and learn from biological data._

LaminDB is an open-source data lake for R&D in biology. It manages indexed **object storage** (local directories, S3, GCP) with a mapped **SQL database** (SQLite, Postgres, and soon, BigQuery).

One cool thing is that you can readily create distributed _LaminDB instances_ at any scale. Get started on your laptop, deploy in the cloud, or work with a mesh of instances for different teams and purposes.

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

```

## Installation

LaminDB is a python package available for Python versions 3.8+.

```shell
pip install lamindb
```

## Import

In your python script, import LaminDB as:

```python
import lamindb as ln
```

## Quick setup

Quick setup on the command line:

- Sign up via `lamin signup <email>`
- Log in via `lamin login <handle>`
- Set up an instance via `lamin init --storage <storage> --schema <schema_modules>`

:::{dropdown} Example code

```shell
lamin signup testuser1@lamin.ai
lamin login testuser1
lamin init --storage ./mydata --schema bionty,wetlab
```

:::

See {doc}`/guide/setup` for more.

## Track & query data

### Track data source & data

::::{tab-set}
:::{tab-item} Within a notebook

```{code-block} python
ln.nb.header()  # data source is created and linked

df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})

# create a data object with SQL metadata record
dobject = ln.DObject(df, name="My dataframe")

# upload the data file to the configured storage
# and commit a DObject record to the SQL database
ln.add(dobject)
```

:::
:::{tab-item} Within a pipeline

```{code-block} python
# create a pipeline record
pipeline = lns.Pipeline(name="my pipeline", version="1")

# create a run from the above pipeline as the data source
run = lns.Run(pipeline=pipeline, name="my run")

df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})

# create a data object with SQL metadata record
dobject = ln.DObject(df, name="My dataframe", source=run)

# upload the data file to the configured storage
# and commit a DObject record to the SQL database
ln.add(dobject)
```

:::
::::

### Query & load data

```python
dobject = ln.select(ln.DObject, name="My dataframe").one()
df = dobject.load()
```

<br>

See {doc}`/guide/track` for more.

## Track biological features

```python
import bionty as bt

# An sample single cell RNA-seq dataset
adata = ln.dev.datasets.anndata_mouse_sc_lymph_node()

# Start to track genes mapped to a Bionty Entity
# - ensembl id as the standardized id
# - mouse as the species
reference = bt.Gene(id=bt.gene_id.ensembl_gene_id, species=bt.Species().lookup.mouse)

# Create a data object with features
dobject = ln.DObject(adata, name="Mouse Lymph Node scRNA-seq", features_ref=reference)

# upload the data file to the configured storage
# and commit a DObject record to the sql database
ln.add(dobject)
```

<br>

See {doc}`/guide/link-features` for more.

```{tip}
- Each page in this guide is a Jupyter Notebook, which you can download [here](https://github.com/laminlabs/lamindb/tree/main/docs/guide).
- You can run these notebooks in hosted versions of JupyterLab, e.g., [Saturn Cloud](https://github.com/laminlabs/run-lamin-on-saturn), Google Vertex AI, and others.
- We recommend using [JupyterLab](https://jupyterlab.readthedocs.io/) for best notebook tracking experience.
```

📬 [Reach out](https://lamin.ai/contact) to report issues, learn about data modules that connect your assays, pipelines & workflows within our data platform enterprise plan.
