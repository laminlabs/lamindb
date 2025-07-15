[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB - A data framework for biology

LaminDB is an open-source data framework to enable learning at scale in computational biology.
It lets you track data transformations, curate datasets, manage metadata, and query a built-in database for biological entities & data structures.

## Setup

<!-- quick-setup-lamindb.md -->

Install the `lamindb` Python package:

```shell
pip install 'lamindb[jupyter,bionty]'  # support notebooks & biological ontologies
```

Create a local LaminDB instance:

```shell
lamin init --storage ./quickstart-data
```

Or if you have write access to an existing instance, simply connect to it:

```shell
lamin connect account/name
```

## Quickstart

<!-- py-quickstart.py -->

```python
import lamindb as ln

ln.track()  # track the run of a script or notebook
open("sample.fastq", "w").write("@r1\nACGT\n+\nIIII\n")  # create some data
ln.Artifact("sample.fastq", key="sample.fastq").save()  # create a versioned artifact
ln.finish()  # finish the run, save source code & run report
```

## Docs

- Copy [summary.md](https://docs.lamin.ai/summary.md) into an LLM chat and let AI explain LaminDB.
- Read the comprehensive [docs](https://docs.lamin.ai).
