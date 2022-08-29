```{include} ../README.md
:start-line: 0
:end-line: 4
```

```{warning}

Public beta: Support is currently only offered to partners & enterprise customers.

```

Curate, store, track, query, integrate, and learn from biological data.

Modular configurable data & analysis platform for hybrid R&D organizations to

<!-- Key features/functionality/capabilities -->

1. query low- and high-dimensional data by biological entities
2. query data by provenance (users, notebooks, pipelines, instruments, etc.)
3. collaborate on data within and across organizations

<!-- User experience -->

with

1. an intuitive API to connect data and analytics infrastructure
2. zero lock-in danger due to an open-source & multi-cloud stack
3. a tool to easily manage schema module migrations in a changing R&D environment
4. support for learning from data across measured → relevant → derived features
5. support for fast-paced iterations and "development data" through data versioning, quality & integrity flags

<!-- High-level technical specification -->

LaminDB is a distributed data management system similar to how git is a distributed version control system.
Each LaminDB instance is a data warehouse with storage (local directory, S3, GCP, Azure) and a SQL database (SQLite, Postgres, BigQuery) for querying it.

Install:

```
pip install lamindb
```

Get started:

- [Tutorials](tutorials/index) walk you through setup and usage of the platform.
- Explore real-world [examples](examples/index).
- Browse the [API reference](api).
- If you get stuck, see [guides](guides/index) for edge cases & errors.

References:

- See [lamin.ai/docs](https://lamin.ai/docs) for an overview of associated open-source modules.
- [Reach out](https://lamin.ai/contact) to learn about modules that connect your assays, pipelines, instruments & workflows within our data platform enterprise offer.
- Read the following reports to learn about technology underlying LaminDB: ...

```{toctree}
:maxdepth: 1
:hidden:

tutorials/index
examples/index
api
guides/index
notes/index
changelog
```
