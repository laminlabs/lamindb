[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB: Manage R&D data & analyses

LaminDB is an open-source data lake for biology.

It gives you components to build on data lineage & biological entities with an ORM for your existing infrastructure: **object storage** (local directories, S3, GCP) with a mapped **SQL query engine** (currently: SQLite, Postgres).

You can create distributed **LaminDB instances** at any scale:

- Get started on your laptop, deploy in the cloud, or work with a mesh of instances for different teams and purposes.
- Share them through a hub akin to HuggingFace & GitHub - see, e.g, [lamin.ai/sunnyosun](https://lamin.ai/sunnyosun).

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

Update 2023-06-05: We completed a major migration from SQLAlchemy/SQLModel to Django, available in pre-release of v0.42.

```

## Installation

LaminDB is a Python package (Python 3.8+).

```shell
pip install lamindb
```

<br>

It is configurable & works with custom schema modules (each being managed as a Django app):

```shell
pip install 'lamindb[bionty]'  # install biological entities
pip install 'lamindb[nbproject]'  # install Jupyter notebook tracking
pip install 'lamindb[aws]'  # install AWS dependencies (s3fs, etc.)
pip install 'lamindb[gcp]'  # install GCP dependencies (s3fs, etc.)
```

<br>

## Setup

Quick setup on the command line:

- Sign up via `lamin signup <email>`
- Log in via `lamin login <handle>`
- Set up an instance via `lamin init --storage <storage> --schema <schema_modules>`

## Documentation

Read the [docs](https://lamin.ai/docs/guide/).
