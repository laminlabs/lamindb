[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB: Data lake for biology

LaminDB is an API layer for your existing infrastructure to manage your existing data & analyses.

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

Update 2023-06-05: We completed a major migration from SQLAlchemy/SQLModel to Django, available in pre-releases of v0.42.

```

## Features

Free:

- Track data lineage across notebooks, pipelines & apps.
- Manage biological registries, ontologies & features.
- Persist, load & stream data objects with a single line of code.
- Query for anything, define & manage custom schemas.
- Manage data on your laptop, on your server or in your cloud infra.
- Use a mesh of distributed LaminDB instances for different teams and purposes.
- Share instances through a Hub akin to GitHub.

Enterprise:

- Explore & share data, submit samples & track lineage with LaminApp (deployable in your infra).
- Receive support, code templates & services for a BioTech data & analytics platform.

## How does it work?

LaminDB builds semantics of R&D and biology onto well-established tools:

- SQLite & Postgres for SQL databases using Django ORM (previously: SQLModel)
- S3, GCP & local storage for object storage using fsspec
- Configurable storage formats: pyarrow, anndata, zarr, etc.
- Biological knowledge sources & ontologies: see [Bionty](https://lamin.ai/docs/bionty)

LaminDB is open source. For details, see [Architecture](#architecture).

## Installation

![pyversions](https://img.shields.io/pypi/pyversions/lamindb)

```shell
pip install lamindb  # basic data lake
pip install 'lamindb[bionty]'  # biological entities
pip install 'lamindb[nbproject]'  # Jupyter notebook tracking
pip install 'lamindb[aws]'  # AWS dependencies (s3fs, etc.)
pip install 'lamindb[gcp]'  # GCP dependencies (gcfs, etc.)
```

## Quick setup

Why do I have to sign up?

- Data lineage requires a user identity (who modified which data when?).
- Collaboration requires a user identity (who shares this with me?).

Signing up takes 1 min.

We do _not_ store any of your data, but only basic metadata about you (email address, etc.) & your instances (S3 bucket names, etc.).

- Sign up: `lamin signup <email>`
- Log in: `lamin login <handle>`
- Init a data lake: `lamin init --storage <default-storage-location>`

## Usage overview

```
import lamindb as ln
```

### Store & load data artifacts

Store a `DataFrame` or an `AnnData` in the default local or cloud storage location:

```python
df = pd.DataFrame({"feat1": [1, 2], "feat2": [3, 4]})

ln.File(df, name="My dataframe").save()  # create a File object and save it
```

<br>

Get it back:

```python
file = ln.File.select(name="My dataframe").one()  # query for it
df = file.load()  # load it into memory
  feat1 feat2
0     1     3
1     2     4
```

### Track & query data lineage

#### Basics

```python
ln.File.select(created_by__handle="user1").df()   # a DataFrame of all files ingested by user1
ln.File.select().order_by("-updated_at").first()  # latest updated file
```

#### Notebooks

Track a Jupyter Notebook:

```python
ln.track()  # auto-detect a notebook, save as a Transform, create a Run
ln.File("my_artifact.parquet").save()  # this file is an output of the notebook run
```

<br>

When you query this file later on, you'll always know where it came from:

```{python}

file = ln.File.select(name="my_artifact.parquet").one()
file.transform  # gives you the notebook with title, filename, version, id, etc.
file.run  # gives you the run of the notebook that created the file
```

<br>

Of course, you can also query for notebooks & the artifacts they created:

```python
transforms = ln.Transform.select(  # all notebooks with 'T cell' in the title created in 2022
    name__contains="T cell", type="notebook", created_at__year=2022
).all()
ln.File.select(transform__in=transforms).all()  # data artifacts created by these notebooks
```

#### Pipelines

To save a pipeline (complementary to workflow tools) to the `Transform` registry, call

```
ln.Transform(name="Awesom-O", version="0.41.2").save()  # save a pipeline
```

<br>

To track a run of a registered pipeline:

```python
transform = ln.Transform.select(name="Awesom-O", version="0.41.2").one()  # select a pipeline from the registry
ln.track(transform)  # create a new global run context
ln.File("s3://my_samples01/my_artifact.fastq.gz").save()  # link file against run & transform
```

<br>

Now, you can query, e.g., for

```python
ln.select(ln.Run, transform__name="Awesom-O").order_by("-created_at").df()  # get the latest Cell Ranger pipeline runs
```

### Manage biological registries

```
lamin init --storage ./myobjects --schema bionty
```

<br>

...

### Track biological features

...

### Track biological samples

...

### Manage custom schemas

1. Create a GitHub repository with Django ORMs similar to [github.com/laminlabs/lnschema-lamin1](https://github.com/laminlabs/lnschema-lamin1)
2. Create & deploy migrations via `lamin migrate create` and `lamin migrate deploy`

It's fastest if we do this for you based on our templates within an enterprise plan, but you can fully manage the process yourself.

## Notebooks

- Find all guide notebooks [here](https://github.com/laminlabs/lamindb/tree/main/docs/guide).
- You can run these notebooks in hosted versions of JupyterLab, e.g., [Saturn Cloud](https://github.com/laminlabs/run-lamin-on-saturn), Google Vertex AI, and others or on Google Colab.
- Jupyter Lab & Notebook offer a fully interactive experience, VS Code & others require using the CLI (`lamin track my-notebook.ipynb`)

## Architecture

LaminDB consists of the `lamindb` Python package, which builds on a number of open-source packages:

- [bionty](https://github.com/laminlabs/bionty): Basic biological entities (usable standalone).
- [lamindb-setup](https://github.com/laminlabs/lamindb-setup): Setup & configure LaminDB, client for Lamin Hub.
- [lnschema-core](https://github.com/laminlabs/lnschema-core): Core schema, ORMs to model data objects & data lineage.
- [lnschema-bionty](https://github.com/laminlabs/lnschema-bionty): Bionty schema, ORMs that are coupled to Bionty's entities.
- [lnschema-lamin1](https://github.com/laminlabs/lnschema-lamin1): Exemplary configured schema to track samples, treatments, etc.
- [nbproject](https://github.com/laminlabs/nbproject): Parse metadata from Jupyter notebooks.

LaminHub & LaminApp are not open-sourced, neither are templates that model lab operations.

Lamin's packages build on the infrastructure listed
[above](#how-does-it-work). Previously, they were based on SQLAlchemy/SQLModel
instead of Django, and cloudpathlib instead of fsspec.

## Documentation

Read the [docs](https://lamin.ai/docs/guide/).
