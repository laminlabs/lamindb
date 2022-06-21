```{include} ../README.md
:start-line: 0
:end-line: 1
```

Interactively curate, store, track and query biological data.

Start out with unstructured collections of files and derive a robust data warehouse, built for queries and learning.[^files]

Leverage both a scientific (biological entities) and a technical (data provenance) framework to organize data.

[^files]: Biology can't get rid of files entirely: biological systems are characterized by measurement devices that generate files with data.

## Features

- Query & standardize against biological entities ([Bionty](https://lamin.ai/bionty)) to avoid semantic mess & data silos through ad hoc conventions.
- Full data provenance on the pipeline **and** data science level using [nbproject](https://lamin.ai/nbproject).
- Balance pragmatic, fast-paced R&D operations and robustness. Allow partial annotations, and flag data with an ontology of quality & integrity flags.
- Leverage a delightful default schema that you can customize to your needs!
- Represent the learning process in the schema and how you think about querying data: `measured_features` -> `relevant_features` -> `derived_features`.
- Monitor a full access log to the database: all updates are tracked, nothing ever gets lost again!

## Get started

Build your learning & analytics workflows and any application on top!

- Get started with the [tutorial](tutorial/index).
- Explore and play with real-world [examples](examples/index).
- Browse the full [API reference](api).
- Look up guides [guides](guides/index) that solve specific problems or illustrate common errors.

## Installation and configuration

Install the development version:

```
pip install -e ".[dev,test]"
```

```{toctree}
:maxdepth: 1
:hidden:

tutorial/index
guides/index
examples/index
api
changelog
```
