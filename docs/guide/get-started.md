# Get started

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

```

## What is LaminDB?

LaminDB is a distributed data management system in which users collaborate on DB _instances_.

This is analogous to how developers collaborate on code in repositories, but unlike git and dvc, LaminDB is queryable by entities.[^integrate]

[^integrate]: Like git helps with integrating code across repositories, LaminDB helps with integrating data across instances.

## What is a LaminDB instance?

Each LaminDB instance is

- a [data lakehouse](https://www.databricks.com/glossary/data-lakehouse) that manages indexed object storage (local directories, S3, GCP)
- a SQL query engine (SQLite, Postgres, and soon, BigQuery)

## Features

LaminDB comes with

- data lineage and edit history
- knowledge-managed biological entities for typing and lookups (via [Bionty](https://lamin.ai/docs/bionty))
- configurable R&D [schema modules](https://lamin.ai/docs/db/lamindb.schema)
- tracking of interactive notebooks (via [nbproject](https://lamin.ai/docs/nbproject))

LaminDB is built on open source [Python packages](https://lamin.ai/docs).

```{toctree}
:maxdepth: 1
:hidden:

install
signup
init
schema
```
