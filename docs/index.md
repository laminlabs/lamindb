```{include} ../README.md
:start-line: 0
:end-line: 4
```

```{warning}

Public beta: Support is currently only offered to partners & enterprise customers.

```

Curate, store, track, query, integrate, and learn from biological data.

Data & analysis platform for R&D to

<!-- Key features/functionality/capabilities -->

1. query low- and high-dimensional data by biological entities
2. query data by provenance (users, notebooks, pipelines, instruments, etc.)
3. collaborate on data within and across organizations

<!-- User experience -->

with

1. a unified API across storage and databases
2. zero lock-in due to open-source & multi-cloud
3. schema module migrations
4. configuration of underlying data modules

<!-- High-level technical specification -->

LaminDB is a distributed data management system similar to how git is a distributed version control system.
Each LaminDB instance is a data warehouse with storage (local directory, S3, GCP, Azure) and a SQL database (SQLite, Postgres) for querying it.

Install:

```
pip install lamindb
```

Get started:

- [Tutorials](tutorials/index) walk you through setup and usage of the platform.
- Browse the [API reference](api).
- If you get stuck, see [guides](guides/index) for edge cases & errors.

References:

- See [lamin.ai/docs](https://lamin.ai/docs) for an overview of all open-sourced data modules LaminDB builds on.
- [Reach out](https://lamin.ai/contact) to learn about modules that connect your assays, pipelines, instruments & workflows within our data platform enterprise offer.

```{toctree}
:maxdepth: 1
:hidden:

tutorials/index
api
guides/index
changelog
```
