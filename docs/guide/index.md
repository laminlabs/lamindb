```{include} ../../README.md
:start-line: 0
:end-line: 3
```

# Guide

Welcome to the LaminDB guide! 👋

```{include} ../../README.md
:start-line: 6
:end-line: -4
```

:::{dropdown} Example code

```shell
lamin signup testuser1@lamin.ai
lamin login testuser1
lamin init --storage ./mydata --schema bionty,lamin1
```

:::

## Quickstart

### Track files & metadata with sources

```python
transform = ln.Transform(name="My pipeline", type="pipeline")
ln.track(transform=transform)
# in a notebook, `ln.track()` parses metadata & creates a `Transform` record

df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})

file = ln.File(df, name="My dataframe")
ln.add(file)
```

<br>

Under-the-hood, this created 3 records:

```
Transform(id='OdlFhFWW7qg3', version='0', name='03-memory', type=notebook, title='Track data objects', created_by_id='DzTjkKse', created_at=datetime.datetime(2023, 4, 28, 6, 7, 30))
Run(id='g1xwllJfFZuh24AWKySc', transform_id='OdlFhFWW7qg3', transform_version='0', created_by_id='DzTjkKse', created_at=datetime.datetime(2023, 4, 28, 6, 7, 30))
File(id='DY9JINrVH6sMtqEirMpM', name='iris', suffix='.parquet', size=5629, hash='jUTdERuqlGv_GyqFfIEb2Q', run_id='g1xwllJfFZuh24AWKySc', transform_id='OdlFhFWW7qg3', transform_version='0', storage_id='GLWpJhvg', created_at=datetime.datetime(2023, 4, 28, 6, 7, 32), created_by_id='DzTjkKse')
```

### Query & load files

```python
file = ln.select(ln.File, name="My dataframe").one()
df = file.load()
#>      a	b
#>  0	1	3
#>  1	2	4
```

<br>

Get the file ingested by the latest run:

```python
run = ln.select(ln.Run).order_by(ln.Run.created_at.desc()).first()
file = ln.select(ln.File, run_id=run.id).all()
```

```{tip}
- Each page in this guide is a Jupyter Notebook, which you can download [here](https://github.com/laminlabs/lamindb/tree/main/docs/guide).
- You can run these notebooks in hosted versions of JupyterLab, e.g., [Saturn Cloud](https://github.com/laminlabs/run-lamin-on-saturn), Google Vertex AI, and others.
- We recommend using [JupyterLab](https://jupyterlab.readthedocs.io/) for best notebook tracking experience.
```

📬 [Reach out](https://lamin.ai/contact) to learn about schema modules that connect your assays & workflows within our data platform enterprise plan.

```{toctree}
:hidden:
:caption: Basics

track
select
stream
schema
```

```{toctree}
:hidden:
:caption: Biology

ontologies
features
link-samples
```

```{toctree}
:hidden:
:caption: Other topics

../faq/index
```
