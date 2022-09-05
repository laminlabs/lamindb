```{include} ../README.md
:start-line: 0
:end-line: 4
```

```{warning}

Public beta: Support is currently only offered to partners & enterprise customers.

```

Curate, store, track, query, integrate, and learn from biological data.

Data & analysis platform for R&D to

1. manage low- and high-dimensional data by biological entities
2. manage data by provenance (users, notebooks, pipelines, instruments, etc.)
3. collaborate on data within and across organizations

LaminDB is a distributed data management system similar to how git is a distributed version control system.
Each LaminDB instance is a lakehouse with storage (local directory, S3, GCP, Azure) and a SQL database (SQLite, Postgres) for querying it.

Install:

```
pip install lamindb
```

Get started:

- The [guide](guide/index) walks you through setup and usage of the platform.
- Browse the [API reference](api) and [FAQ](faq/index).

References:

- See [docs](https://lamin.ai/docs) for an overview of open-sourced data modules LaminDB builds on.
- [Reach out](https://lamin.ai/contact) to learn about modules that connect your assays, pipelines & workflows within our data platform enterprise plan.

```{toctree}
:maxdepth: 1
:hidden:

guide/index
api
faq/index
changelog
```
