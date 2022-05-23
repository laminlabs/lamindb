```{include} ../README.md
:start-line: 0
:end-line: 1
```

Interactively & iteratively curate, store, track and query biological data.

Biology can't get rid of files entirely: biological systems will continue to be characterized by measurement devices.
Start out with an unstructured collection of files and end up with a robust data warehouse.

- Designed balancing pragmatic, fast-paced R&D operations and robustness, with data scientist needs in mind.
- Standardize against [Bionty](https://lamin.ai/bionty) to avoid data silos through ad hoc conventions.
- Maintain a mixed in-house & open schema, partial annotations, and annotates data with curation flags.
- Integrate data science workflows using [nbproject](https://lamin.ai/nbproject), to track human-level data provenance.

Build your learning & analytics workflows and any application on top!

```{note}

Currently, the guides & API are in early private beta, heavily re-worked on a weekly basis, and hence yet lack any pedagogics.

```

- Get started with these [guides](guides).
- Browse the full [API reference](api).

## Installation and configuration

Install the development version:

```
pip install -e ".[dev,test]"
```

```{toctree}
:maxdepth: 1
:hidden:

guides
api
```
