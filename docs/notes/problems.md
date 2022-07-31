# Common problems of data-heavy R&D organizations

We'll add the summary here once we completed the detailed table.

## Detailed table

<!-- prettier-ignore -->
| Description | Solution
--- | --- | ---
**Fundamental scientific problems** |
_Stand on solid ground._ | Key analytics results cannot be linked to supporting data as too many processing steps are involved. | LaminDB provides full provenance of analytics results.
_Blob storage._ | Data can't be queried as it's buried in blob storage. | Index observations and variables in blob storage and link them against a query database.
_Pile of data._ | Data can't be accessed as it's not well-structured. | Structure data **both** with the scientific hypothesis space and provenance (data flow and iterative R&D operations).
_Learn from data at scale._ | Data can't be accessed at scale as no viable programmatic interfaces exist. | LaminDB creates an API basis that can power programmatic access and graphical wetlab interfaces.
_Data integration._ | Molecular (high-dimensional) data can't be integrated with phenotypic (low-dimensional) data. | Index molecular data with the same entities as phenotypic data, provide interfaces with databases for low-dimensional data.
**Problems in collaborative science** |
_Data silos: infrastructure._ | Data can't be easily shared across organizations. | Flexible permission system on LaminHub.
_Data silos: semantics._ | External data can't be mapped on in-house data and vice versa. | Provide intuitive curation and ingestion APIs, operate on open-source data models that can be adopted by any organization, LaminHub offers overlapping schema.
**Problems in day-to-day science** |
_Support learning._ | There is no support for the typical learning-from-data cycle | LaminDB supports data models across the full lab cycle, including measured → relevant → derived features, and supports expanding knowledge management through it's rich semantic models along side managing high-dimensional data to learn from
**Problems in day-to-day R&D operations** |
_Development data._ | Development data can't be ingested as data models are too rigid. | LaminDB allows ingesting development data and decorates them with corresponding QC flags to distinguish them from production data.
_Corrupted data._ | Data is corrupted | LaminDB's full provenance allows to trace back corruption to its origin and write a simple fix, typically, in form of an ingestion constraint
**Problems in building data platforms** |
_Aligning on data models._ | Data models are hard to align across interdisciplinary stakeholders. | LaminDB's templates should cover 90% of cases, the remaining 10% can be get configured.
_Lock-in._ | Commercial solutions lock organizations into specific cloud infrastructure | LaminDB works based on an open source and a multi-cloud stack with zero lock-in danger.
_Migrations are a pain._ | Migrating data models in a fast-paced R&D environment is cumbersome. | LaminDB supports schema-module-based migrations.
