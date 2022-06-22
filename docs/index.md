```{include} ../README.md
:start-line: 0
:end-line: 1
```

Interactively curate, store, track, query, and learn from biological data across users and organizations.

Speak SQL & pydata, embrace biological entities & data provenance, and migrate data models & backends with ease.

Features:

- Query by & standardize against [Bionty](https://lamin.ai/bionty) to avoid semantic mess & data silos through ad hoc conventions.
- Query by & track data provenance on the pipeline- & data-science-level using [nbproject](https://lamin.ai/nbproject).
- Balance pragmatic, fast-paced R&D iterations and robustness. Allow partial annotations, flag data with an ontology of quality & integrity flags.
- Easily migrate the default schema to whatever you need. Migrate storage & database backends with ease.
- Ingest any data file from any measurement device and receive guidance for integrating and querying it.
- Built for the alternating cycle of learning from data and querying data: `measured_features` -> `relevant_features` -> `derived_features`.
- Monitor a full access log to the database: everything is tracked, never lose anything again, resolve corrupted data, easily get rid of unused data!
- Integrate data across organizations & users by virtue of universal backends and IDs.

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
