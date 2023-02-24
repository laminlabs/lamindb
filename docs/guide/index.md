# Guide

Welcome to LaminDB! 👋

_Curate, store, track, query, integrate, and learn from biological data._

_LaminDB_ is a distributed data management system in which users collaborate on DB _instances_.

Each _LaminDB instance_ is a [data lakehouse](https://www.databricks.com/glossary/data-lakehouse) that manages indexed **object storage** (local directories, S3, GCP) with a **SQL query engine** (SQLite, Postgres, and soon, BigQuery).

This is analogous to how developers collaborate on code in repositories, but unlike git and dvc, LaminDB is **queryable by entities**.

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

```

## Features

LaminDB comes with

- data lineage and edit history
- tracking of interactive notebooks
- knowledge-managed biological entities for typing and lookups
- configurable schema modules

LaminDB is built on open-source [Python packages](https://lamin.ai/docs).

## Getting started

Quick setup on the command line (see [Initialize a LaminDB instance](https://lamin.ai/docs/db/guide/setup) for advanced setup guide):

- Install via `pip install lamindb`
- Sign up via `lamin signup <email>` and confirm the sign-up email
- Log in via `lamin login <handle>`

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
view
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: Advanced concepts

dobject.md
setup
schema
session
query-book
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: Track runs

run.md
nb
pipeline
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: Track biology

knowledge
link-features
link-samples
select-features
view-bio
```
