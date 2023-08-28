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

View the lineage of a given file (from [here](docs:birds-eye)):

```python
file.view_lineage()
```

<img src="https://raw.githubusercontent.com/laminlabs/lamindb/main/docs/img/readme/view_lineage.svg" width="800">

### Manage biological registries

Create a cell type registry from public knowledge and add a new cell state (from [here](bio-registries)):

```python
import lnschema_bionty as lb  # import basic biological registries

lb.CellType.from_bionty(name="neuron").save()  # create an ontology-coupled cell type record and save it

new_cell = lb.CellType(name="my neuron cell")
neuron = lb.CellType.lookup().neuron  # look up neuron cell type from database
new_cell.parents.add(neuron)  # add neuron as parent
new_cell.view_parents()
```

<img src="https://raw.githubusercontent.com/laminlabs/lamindb/main/docs/img/readme/neuron_view_parents_dist%3D2.svg" width="500">

### Query validated features & labels

Query for rich, validated meta-data and get an overview using `file.describe()` (from [here](docs:scrna)):

```
💡 File(id='nRwgLc3CMIBr84ckw7jj', key=None, suffix='.h5ad', accessor='AnnData', description='Detmar22', version=None, size=17177479, hash='7ni-11ZJqq0h3LKvS2NcyQ', hash_type='md5', created_at=2023-08-24 20:38:51, updated_at=2023-08-24 20:38:51)

Provenance:
    🗃️ storage: Storage(id='FbnPfTQx', root='/home/runner/work/lamin-usecases/lamin-usecases/docs/test-scrna', type='local', updated_at=2023-08-24 20:38:26, created_by_id='DzTjkKse')
    💫 transform: Transform(id='Nv48yAceNSh8z8', name='Validate & register scRNA-seq datasets', short_name='scrna', version='0', type=notebook, updated_at=2023-08-24 20:38:50, created_by_id='DzTjkKse')
    👣 run: Run(id='eVOa2WmrpGmBAR4cvJzZ', run_at=2023-08-24 20:38:28, transform_id='Nv48yAceNSh8z8', created_by_id='DzTjkKse')
    👤 created_by: User(id='DzTjkKse', handle='testuser1', email='testuser1@lamin.ai', name='Test User1', updated_at=2023-08-24 20:38:26)
Features:
  var (X):
    🔗 index (10000, bionty.Gene.id): ['f9vsAIUY3qUV', '5uxXp1P1fnh0', 'bk8oApw4bzwT', 'IYxX3I0M3jTw', '0tY7C8m0z0EM'...]
  external:
    🔗 assay (2, bionty.ExperimentalFactor): ['single-cell RNA sequencing', 'adult']
  obs (metadata):
    🔗 cell_type (1, bionty.CellType): ['endothelial cell']
    🔗 developmental_stage (2, bionty.ExperimentalFactor): ['single-cell RNA sequencing', 'adult']
    🔗 species (1, bionty.Species): ['mouse']
    🔗 tissue (1, bionty.Tissue): ['inguinal lymph node']
    🔗 immunophenotype (2, core.Label): ['CD45 positive', 'CD45 negative']
    🔗 genotype (1, core.Label): ['wild type genotype']
    🔗 age (1, core.Label): ['8 to 10 week']
    🔗 sex (1, core.Label): ['female']
```

### Collaborate across a mesh of instances

If provided with access, others can enjoy validated & queryable data by loading your instance via:

```shell
$ lamin load myhandle/myinstance
```

### Manage custom schemas

1. Create a GitHub repository with registries similar to [github.com/laminlabs/lnschema-lamin1](https://github.com/laminlabs/lnschema-lamin1)
2. Create & deploy migrations via `lamin migrate create` and `lamin migrate deploy`

It's fastest if we do this for you based on our templates within an enterprise plan.

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
tutorial1
```

```{toctree}
:hidden:
:caption: "How to"

query-search
validate
bio-registries
schemas
setup
```

```{toctree}
:hidden:
:caption: Other topics

faq
storage
```
