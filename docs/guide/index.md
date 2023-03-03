# Guide

Welcome to the LaminDB guide! 👋

_Curate, store, track, query, integrate, and learn from biological data._

LaminDB provides distributed data management in which users collaborate on _LaminDB instances_.

Each _LaminDB instance_ is a [data lakehouse](https://www.databricks.com/glossary/data-lakehouse) that manages indexed **object storage** (local directories, S3, GCP) with a mapped **SQL query database** (SQLite, Postgres, and soon, BigQuery).

This is analogous to how developers collaborate on code in repositories, but unlike git and dvc, LaminDB is **queryable by entities**.

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

```

## Installation

LaminDB is a python package available for Python versions 3.8+.

```bash

pip install lamindb
```

## Quick sign-up and setup

Quick setup on the command line (see [Initialize a LaminDB instance](https://lamin.ai/docs/guide/setup) for advanced setup guide):

- Sign up via `lamin signup <email>` and confirm the sign-up email
- Log in via `lamin login <handle>`

## Tracking data with LaminDB

Inside a notebook:

```python
# test-lamin.ipynb
import lamindb as ln

# tracks the notebook run as a data source.
ln.nb.header()

filepath = "./myproject/mypic.png"
# start tracking your file
dobject = ln.DObject(filepath)
ln.add(dobject)
```

With a python script:

```python
# test-lamin.py

import lamindb as ln
import lamindb.schema as lns

# create a run from a pipeline as the data source
pipeline = lns.Pipeline(name="my pipeline", version="1")
run = lns.Run(pipeline=pipeline, name="my run")

filepath = "./myproject/mypic.png"
# start tracking your file
dobject = ln.DObject(filepath)
ln.add(dobject)
```

```{tip}

- Each page in this guide is a Jupyter Notebook, which you can download [here](https://github.com/laminlabs/lamindb/tree/main/docs/guide).
- You can run these notebooks in hosted versions of JupyterLab, e.g., [Saturn Cloud](https://github.com/laminlabs/run-lamin-on-saturn), Google Vertex AI, and others.
- We recommend using [JupyterLab](https://jupyterlab.readthedocs.io/) for best notebook tracking experience.

```

📬 [Reach out](https://lamin.ai/contact) to report issues, learn about data modules that connect your assays, pipelines & workflows within our data platform enterprise plan.

```{toctree}
:maxdepth: 1
:hidden:
:caption: Get started

ingest
select
add-delete
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: More details

dobject.md
setup
schema
query-book
run
session
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: Track biology, features, samples

knowledge
link-features
link-samples
```
