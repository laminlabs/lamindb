[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB: Data lakes for biology

LaminDB is an API layer for your existing infrastructure to manage your existing data & analyses.

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

Update 2023-06-05: We completed a major migration from SQLAlchemy/SQLModel to Django, available in pre-releases of v0.42.

```

## Features

Free:

- track data lineage across notebooks, pipelines & apps (& integrate with workflow tools)
- manage biological registries, ontologies & features
- persist, load & stream data objects with a single line of code
- query for anything & everything
- define & manage your own schemas (assays, instruments, etc.)
- manage data on your laptop, on your server or in your cloud infra
- use a mesh of distributed LaminDB instances for different teams and purposes
- share instances through a hub akin to GitHub

If you want more, [reach out](https://lamin.ai/contact) for an enterprise plan to:

- explore & share data, submit samples & track lineage with LaminApp (deployable in your infra)
- receive support & services for a BioTech data & analytics platform

## How does it work?

LaminDB builds semantics of R&D and biology onto well-established tools:

- SQLite & Postgres for SQL databases
- S3, GCP & local storage for object storage
- Django for an ORM
- configurable storage backends to persist data objects: pyarrow, anndata, zarr, etc.
- biological knowledge resources & ontologies: see [Bionty](https://lamin.ai/docs/bionty)

Most of LaminDB is open source.

## Installation

```shell
pip install lamindb  # basic data lake
pip install 'lamindb[bionty]'  # biological entities
pip install 'lamindb[nbproject]'  # Jupyter notebook tracking
pip install 'lamindb[aws]'  # AWS dependencies (s3fs, etc.)
pip install 'lamindb[gcp]'  # GCP dependencies (gcfs, etc.)
```

## Quick setup

Why do I have to sign up?

- Data lineage requires a unique user identity (who modified which data when?)
- Sharing & collaborating on data requires a user identity

Signing up takes 1 min and we only store basic metadata about you (email address, etc.) & your instances (S3 bucket names, etc.).

We don't store any of your data!

- Sign up via `lamin signup <email>`
- Log in via `lamin login <handle>`
- Init an instance via `lamin init --storage <storage>`

## Quick overview

### Track & query data lineage

```python
ln.track()  # auto-detect a notebook & register as a Transform
ln.File("my_artifact.parquet").save()  # link Transform & Run objects to File object
```

These 2 lines enable to

```python
ln.File.select(created_by="user1").df()   # all files ingested by user1
ln.File.select(created_by="user1").order_by("-updated_at").first()   # latest modified file by user1

transforms = ln.Transform.select(name__contains="T cell", type="notebook").all()  # all notebooks with 'T cell' in the title
ln.File.select(transform=transforms[1]).all()  # get the file ingested by the second notebook in transforms
# etc.
```

Or, if you'd like to track a run of a register pipeline (here, "Cell Ranger"):

```python
transform = ln.Transform.select(name="Cell Ranger", version="0.7.1").one()  # select a pipeline from the registry
run = ln.Run(transform)  # create a new run record
ln.File("s3://my_samples01/my_artifact.fastq.gz", run=run).save()  # link file against run
# alternatively, create a global run context ln.track(transform) and omit run=run
```

This enables to:

```python
run = ln.select(ln.Run, transform__name="Cell Ranger").order_by("-created_at").df()  # get the latest Cell Ranger pipeline runs
# query files by selected runs, etc.
```

### Serialize & load data objects

```python
df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})

ln.File(df, name="My dataframe").save()
```

<br>

Get it back:

```python
file = ln.select(ln.File, name="My dataframe").one()  # query for it
df = file.load()  # load it into memory
    a   b
0   1   3
1   2   4
```

### Mange biological registries

```
lamin init --storage ./myobjects --schema bionty
```

...

### Track biological features

...

### Track biological samples

...

### Create a custom schema module & manage migrations

1. Create a GitHub repository with the Django models that match your data mimicking [github.com/laminlabs/lnschema-lamin1](https://github.com/laminlabs/lnschema-lamin1)
2. Create & deploy migrations via `lamin migrate create` and `lamin migrate deploy`

It's fastest if we do this for you based on our templates within an enterprise plan, but you can fully manage the process yourself.

## Notebooks

- Find all guide notebooks [here](https://github.com/laminlabs/lamindb/tree/main/docs/guide).
- You can run these notebooks in hosted versions of JupyterLab, e.g., [Saturn Cloud](https://github.com/laminlabs/run-lamin-on-saturn), Google Vertex AI, and others or on Google Colab.
- Jupyter Lab & Notebook offer a fully interactive experience, VS Code & others require using the CLI (`lamin track my-notebook.ipynb`)

## Architecture

LaminDB consists of the `lamindb` Python package, which builds on a number of open-source packages developed by Lamin:

- [bionty](https://github.com/laminlabs/bionty): Biological entities (usable standalone)
- [lamindb-setup](https://github.com/laminlabs/lamindb-setup): Setup & configure LaminDB, client for Lamin Hub.
- [lnschema-core](https://github.com/laminlabs/lnschema-core): Core schema, containing the core ORMs.
- [lnschema-bionty](https://github.com/laminlabs/lnschema-bionty): Bionty schema, containing ORMs that are coupled to Bionty's biological entities.
- [lnschema-lamin1](https://github.com/laminlabs/lnschema-lamin1): Exemplary configured schema to track samples, treatments, etc.
- [nbproject](https://github.com/laminlabs/nbproject): Parse metadata from Jupyter notebooks.

LaminHub & LaminApp are not open sourced, neither are templates to model lab operations.

## Documentation

Read the [docs](https://lamin.ai/docs/guide/).
