[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB - Open-source data platform for biology

```{warning}

Public beta: Close to having converged a stable API, but some breaking changes might still occur.

Here is an [intro video](https://www.youtube.com/watch?v=DtJ9KnqWA8Q) to guide beta testing.

```

LaminDB is a Python library to manage data & analyses related to biology:

- Query, validate & link data batches using biological registries & ontologies.
- Track data flow across pipelines, notebooks & app uploads.
- Manage features & labels schema-less or schema-full.
- Collaborate across a mesh of LaminDB instances.

If you want a UI: LaminApp is built on LaminDB. If LaminDB ~ git, LaminApp ~ GitHub.

(Enterprise features for LaminApp, support, integration tests & schemas are available on a paid plan - in your or our infrastructure.)

## Quickstart

Run `pip install 'lamindb[jupyter]'` and `lamin signup <email>` on the command line (more [info](https://lamin.ai/docs/setup)).

Init a LaminDB instance with local or cloud default storage like you'd init a git repository:

```shell
$ lamin init --storage ./mydata   # or s3://my-bucket, gs://my-bucket
```

Validate & register a `DataFrame` that comes with basic metadata:

```python
import lamindb as ln
import pandas as pd

ln.track()  # track data flow (here: notebook run context)

# define validation criteria by populating the Feature registry
names_types = [("CD14", int), ("CD45", int), ("perturbation", "category")]
features = [ln.Feature(name=name, type=type) for (name, type) in names_types]
ln.save(features)

# access a new batch of data
df = pd.DataFrame(
    {"CD14": [1, 2, 3], "CD45": [3, 4, 5], "perturbation": ["DMSO", "IFNG", "DMSO"]}
)

# validate features & register a Dataset
dataset = ln.Dataset.from_df(df, name="Immune phenotyping 1")
dataset.save()  # save/upload dataset
```

Search, query, and load a dataset:

```python
ln.Dataset.search("immune")  # run a search

# run a query (under the hood, you have the full power of SQL to query)
dataset = ln.Dataset.filter(name__contains="phenotyping 1").one()

# view data flow that generated the dataset
dataset.view_flow()

# load the dataset
df = dataset.load()
```

Use the `bionty` plug-in to type biological entities. For instance, register a whole panel of cell markers:

```python
import lnschema_bionty as lb

# populate the Cell Marker registry using a public ontology
cell_markers = lb.CellMarker.from_values(["CD14", "CD45"])
ln.save(cell_markers)
# define feature validation criteria by registering a panel of markers
marker_panel = ln.FeatureSet(cell_markers, type=int)
ln.save(marker_panel)

# validate a new batch of data against the marker_panel
dataset = ln.Dataset.from_df(df, name="Immune phenotyping 1", columns_ref=lb.CellMarker.name)
dataset.save()  # register dataset
```

## Documentation

Read the [docs](https://lamin.ai/docs).
