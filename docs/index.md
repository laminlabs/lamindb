```{include} ../README.md
:start-line: 0
:end-line: 4
```

_For background, see [Lamin Blog #3 (2022)](https://lamin.ai/notes/2022/lamindb)._

Interactively curate, store, track, query, and learn from biological data across users and organizations.

Speaks SQL & pydata in queries for biological entities & data provenance. Built for migrations and zero-lock-in.

Features:

- Query & track data provenance for pipeline (TODO) & data-science ([nbproject](https://lamin.ai/nbproject)) workflows.
- Query by canonical biological entities ([Bionty](https://lamin.ai/bionty)) to avoid semantic silos.
- Track all user modifications of data.
- Integrate data across users & organizations through universal IDs, backends and overlapping schema.
- Easily migrate the default schema and storage & database backends.
- Ingest data from any instrument and receive guidance for integrating and querying it (TODO).
- Built for the alternating cycle of learning from data and querying data across measured → relevant → derived features.
- Balance pragmatic, fast-paced R&D iterations and robustness by allowing partial annotations combined with quality & integrity flags.

Install:

```
pip install lamindb
```

Get started:

- [Tutorials](tutorials/index) walks you through setup and usage of the platform.
- Explore and play with real-world [examples](examples/index).
- Browse the full [API reference](api).
- Look up guides [guides](guides/index) that solve specific problems or illustrate common errors.

```{toctree}
:maxdepth: 1
:hidden:

tutorials/index
examples/index
api
guides/index
changelog
```
