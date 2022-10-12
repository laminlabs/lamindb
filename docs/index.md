```{include} ../README.md
:start-line: 0
:end-line: 4
```

_Curate, store, track, query, integrate, and learn from biological data._

```{warning}

Public beta: Support is currently only offered to enterprise customers & collaborators.

```

LaminDB is a distributed data management system similar to how git is a distributed version control system.
Each LaminDB instance is a [lakehouse](https://www.databricks.com/glossary/data-lakehouse) with storage locations (local directories, S3, GCP) and a SQL database (SQLite, Postgres) for selecting them.
LaminDB speaks biology via [Bionty](https://lamin.ai/docs/bionty) and configurable R&D [schema modules](https://lamin.ai/docs/db/lamindb.schema).

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
