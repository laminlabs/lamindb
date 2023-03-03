[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB: Manage R&D data & analyses

_Curate, store, track, query, integrate, and learn from biological data._

**Public beta:** Currently only recommended for collaborators as we still make breaking changes.

Read the **[docs](https://lamin.ai/docs)**.

## Tracking data with LaminDB

Inside a notebook:

```python
# test-lamin.ipynb
import lamindb as ln

# tracks the notebook run as a data source.
ln.nb.header()

# track a local file
filepath = "./myproject/myimage.png"
dobject = ln.DObject(filepath)
ln.add(dobject)
```

With a python script:

```python
# test-lamin.py
import lamindb as ln

# create a run from a pipeline as the data source
pipeline = ln.schema.Pipeline(name="My test pipeline")
run = ln.schema.Run(pipeline=pipeline, name="My test run")

# track a local file
filepath = "./myproject/myimage.png"
dobject = ln.DObject(filepath, source=run)
ln.add(dobject)
```
