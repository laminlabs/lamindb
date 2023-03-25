[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB: Manage R&D data & analyses

_Curate, store, track, query, integrate, and learn from biological data._

LaminDB is an open-source data lake for R&D in biology. It manages indexed **object storage** (local directories, S3, GCP) with a mapped **SQL database** (SQLite, Postgres, and soon, BigQuery).

One cool thing is that you can readily create distributed _LaminDB instances_ at any scale. Get started on your laptop, deploy in the cloud, or work with a mesh of instances for different teams and purposes.

```{warning}

Public beta: Currently only recommended for collaborators as we still make breaking changes.

```

## Installation

LaminDB is a python package available for Python versions 3.8+.

```shell
pip install lamindb
```

<br>

If you need to work with bionty (feature parsing) and wetlab schemas:

```shell
pip install 'lamindb[bionty,wetlab]'
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

:::{dropdown} Example code

```shell
lamin signup testuser1@lamin.ai
lamin login testuser1
lamin init --storage ./mydata --schema bionty,wetlab
```

:::

See {doc}`/guide/setup` for more.

## Track & query data

### Track data sources, data, and metadata

::::{tab-set}
:::{tab-item} Within an interactive notebook

```{code-block} python
import lamindb as ln

ln.Run() # data source (a run record) is created
#> ℹ️ Instance: testuser1/mydata
#> ℹ️ User: testuser1
#> ℹ️ Loaded notebook: Notebook(id='OdlFhFWW7qg3', v='0', name='04-memory', title='Track in-memory data objects', created_by='DzTjkKse', created_at=datetime.datetime(2023, 3, 15, 16, 14, 42))
#> ℹ️ Loaded run:
#> Run(id='L1oBMKW60ndt5YtjRqav', notebook_id='sePTpDsGJRq3', notebook_v='0', created_by='bKeW4T6E', created_at=datetime.datetime(2023, 3, 14, 21, 49, 36))

df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})

# create a data object with SQL metadata record including hash
# link run record
file = ln.File(df, name="My dataframe")
#> File(id='dZvGD7YUKCKG4X4aLd5K', name='My dataframe', suffix='.parquet', size=2240, hash='R2_kKlH1nBGesMdyulMYkA', source_id='L1oBMKW60ndt5YtjRqav', storage_id='wor0ul6c')

# upload serialized version to the configured storage
# commit a File record to the SQL database
ln.add(file)
#> File(id='dZvGD7YUKCKG4X4aLd5K', name='My dataframe', suffix='.parquet', size=2240, hash='R2_kKlH1nBGesMdyulMYkA', source_id='L1oBMKW60ndt5YtjRqav', storage_id='wor0ul6c', created_at=datetime.datetime(2023, 3, 14, 21, 49, 46))
```

:::
:::{tab-item} Within a regular pipeline

```{code-block} python
# create (or query) a pipeline record
pipeline = lns.Pipeline(name="My pipeline")
#> Pipeline(id='fhn5Zydf', v='1', name='My pipeline', created_by='bKeW4T6E')

# create a run from the above pipeline as the data source
run = ln.Run(pipeline=pipeline)
#> Run(id='2aaKWH8dwBE6hnj3n9K9', pipeline_id='fhn5Zydf', pipeline_v='1', created_by='bKeW4T6E')

# access pipeline from run via
print(run.pipeline)
#> Pipeline(id='fhn5Zydf', v='1', name='My pipeline', created_by='bKeW4T6E')

df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})

# create a data object with SQL metadata record including hash and link run record
file = ln.File(df, name="My dataframe", source=run)
#> File(id='dZvGD7YUKCKG4X4aLd5K', name='My dataframe', suffix='.parquet', size=2240, hash='R2_kKlH1nBGesMdyulMYkA', source_id='L1oBMKW60ndt5YtjRqav', storage_id='wor0ul6c')

# Tip: If you work with a single thread, you can pass `global_context=True` to ln.Run(), allowing you to omit source=run

# upload serialized version to the configured storage
# commit a File record to the SQL database
ln.add(file)
#> File(id='dZvGD7YUKCKG4X4aLd5K', name='My dataframe', suffix='.parquet', size=2240, hash='R2_kKlH1nBGesMdyulMYkA', source_id='L1oBMKW60ndt5YtjRqav', storage_id='wor0ul6c', created_at=datetime.datetime(2023, 3, 14, 21, 49, 46))
```

:::
::::

### Query & load data

```python
file = ln.select(ln.File, name="My dataframe").one()
#> [File(id='dZvGD7YUKCKG4X4aLd5K', name='My dataframe', suffix='.parquet', size=2240, hash='R2_kKlH1nBGesMdyulMYkA', source_id='L1oBMKW60ndt5YtjRqav', storage_id='wor0ul6c', created_at=datetime.datetime(2023, 3, 14, 21, 49, 46))]
df = file.load()
#>      a	b
#>  0	1	3
#>  1	2	4
```

Get the data ingested by the latest run:

```python
run = ln.select(ln.Run).order_by(ln.Run.created_at.desc()).first()
#> Run(id='L1oBMKW60ndt5YtjRqav', notebook_id='sePTpDsGJRq3', notebook_v='0', created_by='bKeW4T6E', created_at=datetime.datetime(2023, 3, 14, 21, 49, 36))
file = ln.select(ln.File).where(ln.File.source == run).all()
#> [File(id='dZvGD7YUKCKG4X4aLd5K', name='My dataframe', suffix='.parquet', size=2240, hash='R2_kKlH1nBGesMdyulMYkA', source_id='L1oBMKW60ndt5YtjRqav', storage_id='wor0ul6c', created_at=datetime.datetime(2023, 3, 14, 21, 49, 46))]
```

<br>

See {doc}`/guide/track` for more.

## Track biological metadata

### Track biological features

```python
import bionty as bt  # Lamin's manager for biological knowledge
import lamindb as ln

ln.Run()  # assume we're in a notebook and don't need to pass pipeline_name

# a sample single cell RNA-seq dataset
adata = ln.dev.datasets.anndata_mouse_sc_lymph_node()

# Create a reference
# - ensembl id as the standardized id
# - mouse as the species
reference = bt.Gene(species="mouse")

# parse gene identifiers from data and map on reference
features = ln.Features(adata, reference)
#> 🔶 id column not found, using index as features.
#> ✅ 10000 terms (100.0%) are mapped.
#> 🔶 0 terms (0.0%) are not mapped.
# The result is a hashed feature set record:
print(features)
#> Features(id='2Mv3JtH-ScBVYHilbLaQ', type='gene', created_by='bKeW4T6E')
# genes records can be accessed via:
print(features.genes[:3])
#> [Gene(id='ENSMUSG00000020592', species_id='NCBI_10090'),
#>  Gene(id='ENSMUSG00000034931', species_id='NCBI_10090'),
#>  Gene(id='ENSMUSG00000071005', species_id='NCBI_10090')]

# track data with features
file = ln.File(adata, name="Mouse Lymph Node scRNA-seq", features=features)

# access linked gene references
print(file.features.genes[:3])
#> [Gene(id='ENSMUSG00000020592', species_id='NCBI_10090'),
#>  Gene(id='ENSMUSG00000034931', species_id='NCBI_10090'),
#>  Gene(id='ENSMUSG00000071005', species_id='NCBI_10090')]

# upload serialized data to configured storage
# commit a File record to the SQL database
# commit all linked features to the SQL database
ln.add(file)
#> File(id='VRu0Mg93d5l6NLb4znCD', name='Mouse Lymph Node scRNA-seq', suffix='.h5ad', size=17341245, hash='Qprqj0O23197Ko-VobaZiw', source_id='EB78Sl5KPG6wW6XcOlsm', storage_id='0Xt6BY40', created_at=datetime.datetime(2023, 3, 17, 6, 49, 39))
```

<br>

See {doc}`/guide/features` for more.

```{tip}
- Each page in this guide is a Jupyter Notebook, which you can download [here](https://github.com/laminlabs/lamindb/tree/main/docs/guide).
- You can run these notebooks in hosted versions of JupyterLab, e.g., [Saturn Cloud](https://github.com/laminlabs/run-lamin-on-saturn), Google Vertex AI, and others.
- We recommend using [JupyterLab](https://jupyterlab.readthedocs.io/) for best notebook tracking experience.
```

📬 [Reach out](https://lamin.ai/contact) to report issues, learn about data modules that connect your assays, pipelines & workflows within our data platform enterprise plan.
