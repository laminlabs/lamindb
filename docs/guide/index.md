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

:::{dropdown} Example

```shell
lamin signup testuser1@lamin.ai
lamin login testuser1
lamin init --storage ./mydata --schema bionty,lamin1
```

:::

## Quickstart

### Track files & metadata with sources

```python
ln.track()  # auto-detect a notebook

df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})

file = ln.File(df, name="My dataframe")
ln.save(file)
```

<br>

Under-the-hood, this saved the file to storage & created 3 linked SQL records: `file`, `run`, `transform`.

### Query & load files

```python
file = ln.select(ln.File, name="My dataframe").one()
df = file.load()
    a   b
0   1   3
1   2   4
```

<br>

Get the file ingested by the latest run:

```python
run = ln.select(ln.Run).order_by("-created_at").first()
file = ln.select(ln.File, run=run).all()
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

files-folders
provenance
select
stream
schema
```

```{toctree}
:hidden:
:caption: Biology

../biology/ontologies
../biology/features
../biology/link-samples
```

```{toctree}
:hidden:
:caption: Other topics

../faq/index
../storage/index
```
