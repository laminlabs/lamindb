[![Stars](https://img.shields.io/github/stars/laminlabs/lamindb?logo=GitHub&color=yellow)](https://github.com/laminlabs/lamindb)
[![codecov](https://codecov.io/gh/laminlabs/lamindb/branch/main/graph/badge.svg?token=VKMRJ7OWR3)](https://codecov.io/gh/laminlabs/lamindb)
[![pypi](https://img.shields.io/pypi/v/lamindb?color=blue&label=pypi%20package)](https://pypi.org/project/lamindb)

# LaminDB - A data framework for biology

<-- first two sentences sync from preface.md -->

LaminDB is an open-source data framework to enable learning at scale in computational biology.
It lets you track data transformations, curate datasets, manage metadata, and query a built-in database for biological entities & data structures.

## Docs

Copy [summary.md](https://docs.lamin.ai/summary.md) into an LLM chat and let AI explain LaminDB or read the [docs](https://docs.lamin.ai).

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

Here's how to create an artifact while tracking source code, environment, logs, inputs, and outputs of a script or notebook.

<!-- copied from py-quickstart.py -->

```python
import lamindb as ln

ln.track()  # track a run
open("sample.fasta", "w").write(">seq1\nACGT\n")
ln.Artifact("sample.fasta", key="sample.fasta").save()  # create an artifact
ln.finish()  # finish the run, save source code & run report
```

<!-- from here on, slight deviation from preface.md, where all this is treated in the walk through in more depth -->

Running this code inside a script or notebook, e.g., via `python create-fasta.py`, produces the following data lineage.

```python
artifact = ln.Artifact.get(key="sample.fasta")
artifact.view_lineage()
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/EkQATsQL5wqC95Wj0001.png" width="250">

You'll know how that artifact was created.

```python
artifact.describe()
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/BOTCBgHDAvwglN3U0000.png" width="850">

And you can query the artifact by how it was created.

```python
ln.Artifact.filter(transform__key="create-fasta.py").df()
#>                      uid           key                    hash  run_id
#> id
#> 2   4TUnaqJPIJRdsqg60000  sample.fasta  VPvs-qQxRsFFALP6wOgUbg       1
```

Beyond tracking data lineage, LaminDB enables managing datasets in the context of any metadata.

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/6sofuDVvTANB0f480000.png" width="700">
