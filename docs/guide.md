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

### Track & query data flow

View the flow of a given file or dataset (from [here](docs:project-flow)):

```python
file.view_flow()
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

### Validate & query data & metadata

Query for rich, validated meta-data. Get an overview using `file.describe()` (from [here](docs:scrna)):

```
💡 File(id='qrjrjSL5cIOVK4zSUINn', suffix='.h5ad', accessor='AnnData', description='Conde22', size=28049505, hash='WEFcMZxJNmMiUOFrcSTaig', hash_type='md5', created_at=2023-08-28 21:34:02, updated_at=2023-08-28 21:34:02)

Provenance:
    🗃️ storage: Storage(id='U89WhwzI', root='/home/runner/work/lamin-usecases/lamin-usecases/docs/test-scrna', type='local', updated_at=2023-08-28 21:33:26, created_by_id='DzTjkKse')
    💫 transform: Transform(id='Nv48yAceNSh8z8', name='Validate & register scRNA-seq datasets', short_name='scrna', version='0', type=notebook, updated_at=2023-08-28 21:33:57, created_by_id='DzTjkKse')
    👣 run: Run(id='jxoUKhAoN73coXjQZLcG', run_at=2023-08-28 21:33:28, transform_id='Nv48yAceNSh8z8', created_by_id='DzTjkKse')
    👤 created_by: User(id='DzTjkKse', handle='testuser1', email='testuser1@lamin.ai', name='Test User1', updated_at=2023-08-28 21:33:26)
Features:
  var (X):
    🔗 index (36503, bionty.Gene.id): ['yqq0AoR40Zds', 'Fhp0Xl5AlfIA', '7Z6yFaHjoiCK', 'LfFNaOvjysQC', '8PGr70IVfonQ'...]
  external:
    🔗 species (1, bionty.Species): ['human']
  obs (metadata):
    🔗 cell_type (32, bionty.CellType): ['lymphocyte', 'plasmacytoid dendritic cell', 'regulatory T cell', 'macrophage', 'classical monocyte']
    🔗 assay (4, bionty.ExperimentalFactor): ["10x 3' v3", "10x 5' v2", "10x 5' v1", 'single-cell RNA sequencing']
    🔗 tissue (17, bionty.Tissue): ['ileum', 'sigmoid colon', 'transverse colon', 'jejunal epithelium', 'bone marrow']
    🔗 donor (12, core.Label): ['582C', 'A52', '640C', 'A36', 'D503']
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

- SQLite & Postgres for SQL databases using the Django ORM (previously: SQLModel)
- S3, GCP & local storage for object storage using fsspec
- Configurable storage formats: pyarrow, anndata, zarr, etc.
- Biological knowledge sources & ontologies: see [Bionty](https://lamin.ai/docs/bionty)

LaminDB is open source.

### Architecture

LaminDB consists of the `lamindb` Python package (repository [here](https://github.com/laminlabs/lamindb)) with its components:

- [bionty](https://github.com/laminlabs/bionty): Basic biological entities (usable standalone).
- [lamindb-setup](https://github.com/laminlabs/lamindb-setup): Setup & configure LaminDB, client for Lamin Hub.
- [lnschema-core](https://github.com/laminlabs/lnschema-core): Core schema, ORMs to model data objects & data flow.
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
