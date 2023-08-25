[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB - Open-source data platform for biology

```{warning}

Public beta: Close to having converged a stable API, but some breaking changes might still occur.

```

LaminDB is a Python library to manage data & analyses related to biology:

- Track & query data lineage across pipelines, notebooks & app uploads.
- Query, validate & link data batches using biological registries & ontologies.
- Manage features & labels schema-less or schema-full.
- Collaborate across a mesh of LaminDB instances.

If you want a UI: LaminApp is built on LaminDB. If LaminDB ~ git, LaminApp ~ GitHub.

(LaminApp, support, integration tests & schemas for an enterprise platform are available on a paid plan - in your or our infrastructure.)

## Quickstart

[Installation and sign-up](https://lamin.ai/docs/setup) take no time: Run `pip install lamindb` and `lamin signup <email>` on the command line.

Init a LaminDB instance with local or cloud default storage like you'd init a git repository:

```shell
$ lamin init --storage ./mydata   # or s3://my-bucket, gs://my-bucket
```

Validate & register a `DataFrame` that comes with basic metadata:

```python
import lamindb as ln
import pandas as pd

ln.track()  # track run context in a notebook

# save target feature names in Feature registry
features = ln.Feature.from_values(["feature1", "feature2", "perturbation"])
ln.save(features)

# receive a batch of data
df = pd.DataFrame(
    {"feature1": [1, 2, 3], "feature2": [3, 4, 5], "perturbation": ["pert1", "pert2", "pert1"]}
)

# validate features & create a Dataset object
dataset = ln.Dataset.from_df(df, name="Dataset 1")
dataset.save()  # save/upload dataset
```

Search, query, and load a `DataFrame`:

```python
ln.Dataset.search("dataset 1")  # run a search

# run a query (under the hood, you have the full power of SQL to query)
dataset = ln.Dataset.filter(name__contains="set 1").one()

df = dataset.load()
```

## Documentation

Read the [docs](https://lamin.ai/docs/guide/).
