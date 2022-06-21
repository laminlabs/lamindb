```{include} ../README.md
:start-line: 0
:end-line: 1
```

_Interactively curate, store, track, query, and stream biological data._

From collections of files to a data warehouse, built for queries and learning, around biological entities and data science, across teams and organizations.
Biology can't get rid of files entirely: Biological systems are characterized by measurement devices that generate files with data.

- Query & standardize against biological entities ([Bionty](https://lamin.ai/bionty)) to avoid semantic mess & data silos through ad hoc conventions.
- Full data provenance on the pipeline **and** data science level using [nbproject](https://lamin.ai/nbproject).
- Balance pragmatic, fast-paced R&D operations and robustness. Allow partial annotations, and flag data with an ontology of quality & integrity flags.
- Leverage a delightful default schema that you can customize to your needs!
- Represent the learning process in the schema and how you think about querying data: `measured_features` -> `relevant_features` -> `derived_features`.
- Monitor a full access log to the database: all updates are tracked, nothing ever gets lost again!

Install:

```
pip install lamindb
```

Get started:

- [Quickstart](tutorial/quickstart) walks you through setting up a basic platform.
- Explore and play with real-world [examples](examples/index).
- Browse the full [API reference](api).
- Look up guides [guides](guides/index) that solve specific problems or illustrate common errors.

```{toctree}
:maxdepth: 1
:hidden:

tutorial/index
guides/index
examples/index
api
changelog
```
