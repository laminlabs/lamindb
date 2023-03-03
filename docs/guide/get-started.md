# Get started

## Installation

LaminDB is a python package available for Python versions 3.8+.

```bash

pip install lamindb
```

## Quick sign-up and setup

Quick setup on the command line (see [Initialize a LaminDB instance](https://lamin.ai/docs/guide/setup) for advanced setup guide):

- Sign up via `lamin signup <email>` and confirm the sign-up email
- Log in via `lamin login <handle>`

## Tracking data with LaminDB

Inside a notebook:

```python
# test-lamin.ipynb
import lamindb as ln

# tracks the notebook run as a data source.
ln.nb.header()

filepath = "./myproject/mypic.png"
# start tracking your file
dobject = ln.DObject(filepath)
ln.add(dobject)
```

With a python script:

```python
# test-lamin.py

import lamindb as ln
import lamindb.schema as lns

# create a run from a pipeline as the data source
pipeline = lns.Pipeline(name="my pipeline", version="1")
run = lns.Run(pipeline=pipeline, name="my run")

filepath = "./myproject/mypic.png"
# start tracking your file
dobject = ln.DObject(filepath)
ln.add(dobject)
```

```{tip}

- Each page in this guide is a Jupyter Notebook, which you can download [here](https://github.com/laminlabs/lamindb/tree/main/docs/guide).
- You can run these notebooks in hosted versions of JupyterLab, e.g., [Saturn Cloud](https://github.com/laminlabs/run-lamin-on-saturn), Google Vertex AI, and others.
- We recommend using [JupyterLab](https://jupyterlab.readthedocs.io/) for best notebook tracking experience.

```
