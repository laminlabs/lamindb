[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB: Manage R&D data & analyses

_Curate, store, track, query, integrate, and learn from biological data._

LaminDB is an open-source data lake for R&D in biology.

It gives you components to build on data lineage & biological entities with an ORM for your existing infrastructure: **object storage** (local directories, S3, GCP) with a mapped **SQL query engine** (SQLite, Postgres, and soon, BigQuery).

You can readily create distributed **LaminDB instances** at any scale:

- Get started on your laptop, deploy in the cloud, or work with a mesh of instances for different teams and purposes.
- Share them through a hub akin to HuggingFace & GitHub - see, e.g, [lamin.ai/sunnyosun](https://lamin.ai/sunnyosun).

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

```

## Installation

LaminDB is a python package available for Python versions 3.8+.

```shell
pip install lamindb
```

<br>

Biological entities are installed like so:

```shell
pip install 'lamindb[bionty,lamin1]'
```

## Import

In your python script, import LaminDB as:

```python
import lamindb as ln
```

## Quick setup

Quick setup on the command line:

- Sign up via `lamin signup <email>`
- Log in via `lamin login <handle>`
- Set up an instance via `lamin init --storage <storage> --schema <schema_modules>`

## Documentation

Read the [docs](https://lamin.ai/docs/guide/).
