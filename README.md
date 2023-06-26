[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)
[![Documentation](https://img.shields.io/badge/Documentation-green)](https://lamin.ai/docs/guide/)

# LaminDB

Open-source data lake to manage your existing data in your existing infrastructure.

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

Update 2023-06-14:

- We completed a major migration from SQLAlchemy/SQLModel to Django, available in 0.42.0.
- The last version before the migration is 0.41.2.
```

## Features

Free:

- Track [data lineage](https://lamin.ai/docs/guide/data-lineage) across notebooks, pipelines & apps.
- Manage [biological registries, ontologies & features](https://lamin.ai/docs/biology/registries).
- [Query, search & look up anything](https://lamin.ai/docs/guide/select), [manage & migrate custom schemas](https://lamin.ai/docs/setup/migrate).
- [Persist, load](https://lamin.ai/docs/guide/files-records#in-memory-objects) & [stream data objects](https://lamin.ai/docs/guide/stream) with a single line of code.
- [Idempotent](https://lamin.ai/docs/faq/idempotency) and [ACID](https://lamin.ai/docs/faq/acid) operations.
- Use a mesh of LaminDB instances and [share them in a hub](https://lamin.ai/laminlabs) akin to GitHub.

Enterprise:

- Explore, share data & submit samples with LaminApp (deployable in your infrastructure).
- Receive support, code templates & services for a BioTech data & analytics platform.

## Usage overview

Import `lamindb` and initialize a data lake instance with local or cloud default storage:

```python
import lamindb as ln

ln.setup.init(storage="./mydata")  # or s3://my-bucket, gs://my-bucket, etc.
```

### Store, query, search & load data objects

Store a `DataFrame` in default storage:

```python
df = pd.DataFrame({"feat1": [1, 2], "feat2": [3, 4]})  # AnnData works, too

ln.File(df, name="My dataset1").save()  # create a File object and save/upload it
```

You have the full power of SQL to query for metadata, but the simplest query for a file is:

```python
file = ln.File.select(name="My dataset1").one()  # get exactly one result
```

If you don't have specific metadata in mind, run a search:

```python
ln.File.search("dataset1")
```

Once you queried or searched it, load a file back into memory:

```python
df = file.load()
```

Or get a backed accessor to stream its content from the cloud:

```python
backed = file.backed()  # currently works for AnnData, zarr, HDF5, not yet for DataFrame
```

### Store, query & search files

The same API works for any file:

```python
file = ln.File("s3://my-bucket/images/image001.jpg")  # or a local path
file.save()  # register the file
```

Query by `key` (the relative path within your storage):

```python
file.select(key_startswith="images/").df()  # all files in folder "images/" in default storage
```

### Auto-complete categoricals

When you're unsure about spellings, use a lookup object:

```python
users = ln.User.lookup()
ln.File.select(created_by=users.lizlemon)
```

### Track & query data lineage

In addition to basic provenance information (`created_by`, `created_at`,
`created_by`), you can track which notebooks, pipelines & apps
transformed files.

#### Notebooks

Track a Jupyter Notebook:

```python
ln.track()  # auto-detect & save notebook metadata
ln.File("my_artifact.parquet").save()  # this file is now aware that it was saved in this notebook
```

When you query the file, later on, you'll know from which notebook it came:

```python
file = ln.File.select(name="my_artifact.parquet").one()  # query for a file
file.transform  # the notebook with id, title, filename, version, etc.
file.run  # the specific run of the notebook that created the file
```

Alternatively, you can query for notebooks and find the files written by them:

```python
transforms = ln.Transform.select(  # all notebooks with 'T cell' in the title created in 2022
    name__contains="T cell", type="notebook", created_at__year=2022
).all()
ln.File.select(transform__in=transforms).df()  # the files created by these notebooks
```

#### Pipelines

This works like for notebooks just that you need to provide pipeline metadata yourself.

To save a pipeline to the `Transform` registry, call

```python
ln.Transform(name="Awesom-O", version="0.41.2").save()  # save a pipeline, optionally with metadata
```

Track a pipeline run:

```python
transform = ln.Transform.select(name="Awesom-O", version="0.41.2").one()  # select pipeline from the registry
ln.track(transform)  # create a new global run context
ln.File("s3://my_samples01/my_artifact.fastq.gz").save()  # file gets auto-linked against run & transform
```

Now, you can query for the latest pipeline runs:

```python
ln.Run.select(transform=transform).order_by("-created_at").df()  # get the latest pipeline runs
```

#### Run inputs

To track run inputs, pass `is_run_input` to any `File` accessor: `.stage()`, `.load()` or `.backed()`. For instance,

```python
file.load(is_run_input=True)
```

You can also track inputs by default by setting `ln.settings.track_run_inputs = True`.

### Load your data lake from anywhere

If provided with access, others can load your data lake via a single line:

```
$ lamin load myaccount/myartifacts
```

### Manage biological registries

```shell
lamin init --storage ./bioartifacts --schema bionty
```

...

### Track biological features

...

### Track biological samples

...

### Manage custom schemas

1. Create a GitHub repository with Django ORMs similar to [github.com/laminlabs/lnschema-lamin1](https://github.com/laminlabs/lnschema-lamin1)
2. Create & deploy migrations via `lamin migrate create` and `lamin migrate deploy`

It's fastest if we do this for you based on our templates within an enterprise plan, but you can fully manage the process yourself.

## Installation

![pyversions](https://img.shields.io/pypi/pyversions/lamindb)

```shell
pip install lamindb  # basic data lake
pip install 'lamindb[jupyter]'  # Jupyter notebook tracking
pip install 'lamindb[bionty]'  # basic biological entities
pip install 'lamindb[fcs]'  # .fcs files (flow cytometry)
pip install 'lamindb[zarr]'  # zarr storage (streaming arrays)
pip install 'lamindb[aws]'  # AWS (s3fs, etc.)
pip install 'lamindb[gcp]'  # Google Cloud (gcfs, etc.)
```

## Sign up

Why do I have to sign up?

- Data lineage requires a user identity (who modified which data when?).
- Collaboration requires a user identity (who shares this with me?).

Signing up takes 1 min.

We do _not_ store any of your data, but only basic metadata about you (email address, etc.) & your LaminDB instances (S3 bucket names, etc.).

- Sign up: `lamin signup <email>`
- Log in: `lamin login <handle>`

## How does it work?

LaminDB builds semantics of R&D and biology onto well-established tools:

- SQLite & Postgres for SQL databases using Django ORM (previously: SQLModel)
- S3, GCP & local storage for object storage using fsspec
- Configurable storage formats: pyarrow, anndata, zarr, etc.
- Biological knowledge sources & ontologies: see [Bionty](https://lamin.ai/docs/bionty)

LaminDB is open source.

## Architecture

LaminDB consists of the `lamindb` Python package (repository [here](https://github.com/laminlabs/lamindb)) with its components:

- [bionty](https://github.com/laminlabs/bionty): Basic biological entities (usable standalone).
- [lamindb-setup](https://github.com/laminlabs/lamindb-setup): Setup & configure LaminDB, client for Lamin Hub.
- [lnschema-core](https://github.com/laminlabs/lnschema-core): Core schema, ORMs to model data objects & data lineage.
- [lnschema-bionty](https://github.com/laminlabs/lnschema-bionty): Bionty schema, ORMs that are coupled to Bionty's entities.
- [lnschema-lamin1](https://github.com/laminlabs/lnschema-lamin1): Exemplary configured schema to track samples, treatments, etc.
- [nbproject](https://github.com/laminlabs/nbproject): Parse metadata from Jupyter notebooks.

LaminHub & LaminApp are not open-sourced, and neither are templates that model lab operations.

Lamin's packages build on the infrastructure listed [above](#how-does-it-work).

## Notebooks

- Find all guide notebooks [here](https://github.com/laminlabs/lamindb/tree/main/docs/guide).
- You can run these notebooks in hosted versions of JupyterLab, e.g., [Saturn Cloud](https://github.com/laminlabs/run-lamin-on-saturn), Google Vertex AI, Google Colab, and others.
- Jupyter Lab & Notebook offer a fully interactive experience, VS Code & others require using the CLI (`lamin track my-notebook.ipynb`)

## Documentation

Read the [docs](https://lamin.ai/docs/guide/).
