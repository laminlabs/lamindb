[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![Docs](https://img.shields.io/badge/docs-humans-yellow)](https://docs.lamin.ai)
[![DocsLLMs](https://img.shields.io/badge/docs-LLMs-yellow)](https://docs.lamin.ai/summary.md)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)
[![PyPI Downloads](https://img.shields.io/pepy/dt/lamindb?logo=pypi)](https://pepy.tech/project/lamindb)

# LaminDB - A data framework for biology

<!-- first two sentences sync from preface.md -->

LaminDB is an open-source data framework to enable learning at scale in computational biology.
It lets you track data transformations, validate & annotate datasets, and query a built-in database for biological metadata & data structures.

## Setup

<!-- copied from quick-setup-lamindb.md -->

Install the `lamindb` Python package:

```shell
pip install 'lamindb[jupyter,bionty]'  # support notebooks & biological ontologies
```

Create a LaminDB instance:

```shell
lamin init --storage ./quickstart-data  # or s3://my-bucket, gs://my-bucket
```

Or if you have write access to an instance, connect to it:

```shell
lamin connect account/name
```

## Quickstart

<!-- copied from preface.md -->

Track a script or notebook run with source code, inputs, outputs, logs, and environment.

<!-- copied from py-quickstart.py -->

```python
import lamindb as ln

ln.track()  # track a run
open("sample.fasta", "w").write(">seq1\nACGT\n")
ln.Artifact("sample.fasta", key="sample.fasta").save()  # create an artifact
ln.finish()  # finish the run
```

<!-- from here on, slight deviation from preface.md, where all this is treated in the walk through in more depth -->

Running this code inside a script via `python create-fasta.py` produces the following data lineage.

```python
artifact = ln.Artifact.get(key="sample.fasta")  # query artifact by key
artifact.view_lineage()
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/EkQATsQL5wqC95Wj0003.png" width="220">

You'll know how that artifact was created.

```python
artifact.describe()
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/BOTCBgHDAvwglN3U0001.png" width="850">

Conversely, you can query artifacts by the script that created them.

```python
ln.Artifact.get(transform__key="create-fasta.py")  # query artifact by transform key
```

Data lineage is just one type of metadata to help analysis and model training through queries, validation, and annotation. Here is a more [comprehensive example](https://lamin.ai/laminlabs/lamindata/artifact/fgKBV8qdSnbIga0i).

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/6sofuDVvTANB0f480001.png" width="850">

## Docs

Copy [summary.md](https://docs.lamin.ai/summary.md) into an LLM chat and let AI explain or read the [docs](https://docs.lamin.ai).
