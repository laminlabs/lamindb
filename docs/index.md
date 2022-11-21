```{include} ../README.md
:start-line: 0
:end-line: 5
```

_Curate, store, track, query, integrate, and learn from biological data._

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

```

LaminDB is a distributed data management system similar to how git is a distributed version control system.
Each LaminDB instance is a [data lakehouse](https://www.databricks.com/glossary/data-lakehouse) that manages indexed object storage (local directories, S3, GCP) with a SQL query engine (SQLite, Postgres, and soon, BigQuery).

- [Get started](guide/get-started) with the [guide](guide/index), [API](api) and [FAQ](faq/index).
- See basic features on the [landing page](https://lamin.ai).
- See the [source code](https://github.com/laminlabs/lamindb) on GitHub.
- Browse the open-sourced [data modules](https://lamin.ai/docs) LaminDB builds on.

[Reach out](https://lamin.ai/contact) to learn about data modules that connect your assays, pipelines & workflows within our data platform enterprise plan.

```{toctree}
:maxdepth: 1
:hidden:

guide/index
api
faq/index
changelog
```
