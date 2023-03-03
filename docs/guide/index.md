# Guide

Welcome to the LaminDB guide! 👋

_Curate, store, track, query, integrate, and learn from biological data._

LaminDB provides distributed data management in which users collaborate on _LaminDB instances_.

Each _LaminDB instance_ is a [data lakehouse](https://www.databricks.com/glossary/data-lakehouse) that manages indexed **object storage** (local directories, S3, GCP) with a mapped **SQL query database** (SQLite, Postgres, and soon, BigQuery).

This is analogous to how developers collaborate on code in repositories, but unlike git and dvc, LaminDB is **queryable by entities**.

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

```

📬 [Reach out](https://lamin.ai/contact) to report issues, learn about data modules that connect your assays, pipelines & workflows within our data platform enterprise plan.

```{toctree}
:maxdepth: 1
:hidden:
:caption: Get started

get-started.md
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
