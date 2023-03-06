[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB: Manage R&D data & analyses

_Curate, store, track, query, integrate, and learn from biological data._

LaminDB provides distributed data management in which users collaborate on _LaminDB instances_.

Each _LaminDB instance_ is a [data lakehouse](https://www.databricks.com/glossary/data-lakehouse) that manages indexed **object storage** (local directories, S3, GCP) with a mapped **SQL query database** (SQLite, Postgres, and soon, BigQuery).

This is analogous to how developers collaborate on code in repositories, but unlike git and dvc, LaminDB is **queryable by entities**.

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
lamin signup testuser1@ln.setup.ai
lamin login testuser1
lamin init --storage ./mydata --schema bionty,wetlab
```

:::

See {doc}`/guide/setup` for more.

## Tracking data via LaminDB

**To start, create a `Run` object**:
::::{tab-set}
:::{tab-item} Inside a notebook

```{code-block} python
ln.nb.header()

# run will be automatically attached to the data
# run = ln.nb.run
```

:::
:::{tab-item} From a pipeline

```{code-block} python
# create a pipeline record
pipeline = lns.Pipeline(name="my pipeline", version="1")

# create a run from the above pipeline as the data source
run = lns.Run(pipeline=pipeline, name="my run")
```

:::
::::
See {doc}`/guide/run` for more.

**Track data on storage**:
::::{tab-set}
:::{tab-item} Inside a notebook

```{code-block} python
---
emphasize-lines: 5
---
# a file in your local storage
filepath = "./myproject/mypic.png"

# create a data object with sql record and storage
dobject = ln.DObject(filepath)

# upload the data file to the configured storage
# and commit a DObject record to the sql database
ln.add(dobject)
```

:::
:::{tab-item} From a pipeline

```{code-block} python
---
emphasize-lines: 5
---
# a file in your local storage
filepath = "./myproject/mypic.png"

# create a data object with sql record and storage
dobject = ln.DObject(filepath, source=run)

# upload the data file to the configured storage
# and commit a DObject record to the sql database
ln.add(dobject)
```

:::
::::
See {doc}`/guide/ingest` for more.

```{tip}

- Each page in this guide is a Jupyter Notebook, which you can download [here](https://github.com/laminlabs/lamindb/tree/main/docs/guide).
- You can run these notebooks in hosted versions of JupyterLab, e.g., [Saturn Cloud](https://github.com/laminlabs/run-lamin-on-saturn), Google Vertex AI, and others.
- We recommend using [JupyterLab](https://jupyterlab.readthedocs.io/) for best notebook tracking experience.

```

📬 [Reach out](https://lamin.ai/contact) to report issues, learn about data modules that connect your assays, pipelines & workflows within our data platform enterprise plan.
