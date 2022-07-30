```{include} ../README.md
:start-line: 0
:end-line: 4
```

Interactively curate, store, track, query, and learn from biological data.

A modular configurable data & analysis platform for hybrid R&D organizations to

<!-- Key features/functionality/capabilities -->

1. query low- and high-dimensional data by biological entities = "organize data in the hypotheses space"
2. query data by provenance (users, notebooks, pipelines, instruments, etc.) = “track it all”
3. share data within and across organizations in an interoperable, reusable way = “no cleaning & curating anymore”

<!-- User experience -->

It comes with:

1. an intuitive ingestion & query API, which can connect with other data and analytics infrastructure
2. zero lock-in danger due to an open-source & multi-cloud stack
3. a tool to easily manage migrations in a rapidly changing R&D environment
4. support for learning from data across measured → relevant → derived features
5. support for fast-paced R&D iterations and "draft data" through data versioning, quality & integrity flags

<!-- High-level technical specification -->

On the highest technical level, LaminDB is a distributed data management system similar to how git is a distributed version control system.
Each LaminDB instance is a data warehouse with storage (local directory, S3, GCP, Azure) and a SQL database (SQLite, Postgres, BigQuery) for querying it.

Install:

```
pip install lamindb
```

Get started:

- [Tutorials](tutorials/index) walks you through setup and usage of the platform.
- Explore and play with real-world [examples](examples/index).
- Browse the full [API reference](api).
- Look up guides [guides](guides/index) that solve specific problems or illustrate common errors.

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
