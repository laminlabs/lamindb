```{include} ../README.md
:start-line: 0
:end-line: 3
```

# Guide

```{include} ../README.md
:start-line: 6
:end-line: -4
```

## Overview

### Track & query data lineage

View all parent transforms and files (from [here](docs:birds-eye)):

```python
file.view_lineage()
```

<img src="https://raw.githubusercontent.com/laminlabs/lamindb/main/docs/img/readme/view_lineage.svg" width="800">

Track a Jupyter Notebook:

```python
ln.track()  # auto-detect notebook metadata, create global run context
```

Or track a pipeline:

```python
ln.Transform(name="Cell Ranger", version="0.7.2").save()  # save a pipeline to the transform registry

transforms = ln.Transform.filter(name="Cell Ranger", version="0.41.2").one()  # query pipeline from the registry
ln.track(transform)  # create a new global run context
```

### Manage biological registries

```shell
$ lamin init --storage ./biodata --schema bionty  # this time, mout plug in lnschema_bionty
```

```python
import lamindb as ln
import lnschema_bionty as lb  # import basic biological registries

lb.CellType.from_bionty(name="neuron").save()  # create an ontology-coupled cell type record and save it

new_cell = lb.CellType(name="my neuron cell")
neuron = lb.CellType.lookup().neuron  # look up neuron cell type from database
new_cell.parents.add(neuron)  # add neuron as parent
new_cell.view_parents()
```

<img src="https://raw.githubusercontent.com/laminlabs/lamindb/main/docs/img/readme/neuron_view_parents_dist%3D2.svg" width="500">

### Track biological features

```python
# track features present in var(X) and obs
adata = ln.dev.datasets.anndata_with_obs()
file = ln.File.from_anndata(
    adata, description="my RNA-seq dataset", var_ref=lb.Gene.ensembl_gene_id
)
file.save()

# view a summary of tracked features
# you have registered two feature sets: 'obs' and 'var'
file.features

# add labels to features
tissues = lb.Tissue.from_values(adata.obs["tissue"], field=lb.Tissue.name)
diseases = lb.Disease.from_values(adata.obs["disease"], field=lb.Disease.name)
file.features.add_labels(tissues + diseases)

# fetch labels of a feature
file.features["obs"].get(name="tissue").df()

# display rich metadata of a file (provenance and features)
file.describe()
```

### Collaborate across a mesh of instances

If provided with access, others can load your instance via:

```shell
$ lamin load myaccount/myinstance
```

### Manage custom schemas

1. Create a GitHub repository with Django ORMs similar to [github.com/laminlabs/lnschema-lamin1](https://github.com/laminlabs/lnschema-lamin1)
2. Create & deploy migrations via `lamin migrate create` and `lamin migrate deploy`

It's fastest if we do this for you based on our templates within an enterprise plan, but you can fully manage the process yourself.

## Setup

### Installation

![pyversions](https://img.shields.io/pypi/pyversions/lamindb)

```shell
pip install lamindb  # basic data management
```

You can configure the installation using `extras`, e.g.,

```shell
pip install 'lamindb[jupyter,bionty,fcs,aws]'
```

Supported `extras` are:

```yaml
jupyter  # Track Jupyter notebooks
bionty   # Manage basic biological entities
fcs      # Manage FCS files (flow cytometry)
zarr     # Store & stream arrays with zarr
aws      # AWS (s3fs, etc.)
gcp      # Google Cloud (gcfs, etc.)
postgres # Postgres server
```

### Sign up

Why do I have to sign up?

- Data lineage requires a user identity (who modified which data when?).
- Collaboration requires a user identity (who shares this with me?).

Signing up takes 1 min.

We do _not_ store any of your data, but only basic metadata about you (email address, etc.) & your LaminDB instances (S3 bucket names, etc.).

- Sign up: `lamin signup <email>`
- Log in: `lamin login <handle>`

## How does it work?

### Dependencies

LaminDB builds semantics of R&D and biology onto well-established tools:

- SQLite & Postgres for SQL databases using Django Registry (previously: SQLModel)
- S3, GCP & local storage for object storage using fsspec
- Configurable storage formats: pyarrow, anndata, zarr, etc.
- Biological knowledge sources & ontologies: see [Bionty](https://lamin.ai/docs/bionty)

LaminDB is open source.

### Architecture

LaminDB consists of the `lamindb` Python package (repository [here](https://github.com/laminlabs/lamindb)) with its components:

- [bionty](https://github.com/laminlabs/bionty): Basic biological entities (usable standalone).
- [lamindb-setup](https://github.com/laminlabs/lamindb-setup): Setup & configure LaminDB, client for Lamin Hub.
- [lnschema-core](https://github.com/laminlabs/lnschema-core): Core schema, ORMs to model data objects & data lineage.
- [lnschema-bionty](https://github.com/laminlabs/lnschema-bionty): Bionty schema, ORMs that are coupled to Bionty's entities.
- [lnschema-lamin1](https://github.com/laminlabs/lnschema-lamin1): Exemplary configured schema to track samples, treatments, etc.
- [nbproject](https://github.com/laminlabs/nbproject): Parse metadata from Jupyter notebooks.
- [lamin-utils](https://github.com/laminlabs/lamin-utils): Utilities for LaminDB and Bionty.
- [readfcs](https://github.com/laminlabs/readfcs): FCS file reader.

LaminHub & LaminApp are not open-sourced, and neither are templates that model lab operations.

## Notebooks

- Find all guide notebooks [here](https://github.com/laminlabs/lamindb/tree/main/docs/guide).
- You can run these notebooks in hosted versions of JupyterLab, e.g., [Saturn Cloud](https://github.com/laminlabs/run-lamin-on-saturn), Google Vertex AI, Google Colab, and others.
- Jupyter Lab & Notebook offer a fully interactive experience, VS Code & others require using the CLI to track notebooks: `lamin track my-notebook.ipynb`

```{toctree}
:hidden:
:caption: Tutorial

tutorial
```

```{toctree}
:hidden:
:caption: "How to"

meta
data
schemas
setup
bio-registries
lnschema-bionty-overview
```

```{toctree}
:hidden:
:caption: Other topics

faq
storage
```
