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

Create a LaminDB instance:

```shell
lamin init --storage ./quickstart-data  # or s3://my-bucket, gs://my-bucket
```

Or if you have write access to an instance, connect to it:

```shell
lamin connect account/name
```

## Quickstart

Here's how to create an artifact while tracking source code, run environment, run logs, and inputs and outputs of a script or notebook.

<!-- py-quickstart.py -->

```python
import lamindb as ln

ln.track()  # track the run of a script or notebook
open("sample.fastq", "w").write("@r1\nACGT\n+\nIIII\n")
ln.Artifact("sample.fastq", key="sample.fastq").save()  # create a versioned artifact
ln.finish()  # finish the run, save source code & run report
```

<!-- from here on, slight deviation from preface.md, where all this is treated in the walk through in more depth -->

Running the code inside a script or notebook, e.g., via `python create-fastq.py`, produces the following data lineage.

```
artifact = ln.Artifact.get(key="sample.fastq")
artifact.view_lineage()
```

<img src="https://lamin-site-assets.s3.amazonaws.com/.lamindb/EkQATsQL5wqC95Wj0000.png" width="250">

This means you'll always know how that artifact was created.

```python
artifact.describe()
#> Artifact .fastq
#> └── General
#>    ├── uid: 4TUnaqJPIJRdsqg60000          hash: VPvs-qQxRsFFALP6wOgUbg
#>    ├── size: 16 B                         space: all
#>    ├── branch: main                       created_at: 2025-07-15 16:06:25
#>    ├── created_by: falexwolf (Alex Wolf)
#>    ├── key: sample.fastq
#>    ├── storage location / path: /Users/falexwolf/repos/lamin-docs/quickstart-data/.lamindb/4TUnaqJPIJRdsqg60000.fastq
#>    └── transform: py-quickstart.py
```

It also means you can query the artifact by the filename of the script or notebook.

```python
ln.Artifact.filter(transform__key="create-fastq.py").df()
#>                      uid           key                    hash  run_id
#> id
#> 2   4TUnaqJPIJRdsqg60000  sample.fastq  VPvs-qQxRsFFALP6wOgUbg       1
```

## Docs

Copy [summary.md](https://docs.lamin.ai/summary.md) into an LLM chat and let AI explain LaminDB or read the [docs](https://docs.lamin.ai).
