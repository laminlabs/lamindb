# Common problems of data-heavy R&D organizations

## Fundamental scientific problems

<!-- prettier-ignore -->
 Problem | Description | Solution
--- | --- | ---
_Stand on solid ground._ | Key analytics results cannot be linked to supporting data as too many processing steps are involved. | LaminDB provides full provenance of analytics results.
_Blob storage._ | Data can't be queried as it's buried in blob storage. | Index observations and variables in blob storage and link them against a query database.
_Pile of data._ | Data can't be accessed as it's not well-structured. | Structure data both via the scientific hypothesis space and via provenance (data flow and iterative R&D operations).
_Learn from data at scale._ | Data can't be accessed at scale as no viable programmatic interfaces exist. | LaminDB creates an API basis that can power programmatic access and graphical wetlab interfaces.
_Data integration._ | Molecular (high-dimensional) data can't be integrated with phenotypic (low-dimensional) data. | Index molecular data with the same entities as phenotypic data. Provide connectors with storage of low-dimensional data (ELN & LIMS systems).

## Collaborative science

<!-- prettier-ignore -->
 Problem | Description | Solution
--- | --- | ---
_Siloed infrastructure._ | Data can't be easily shared across organizations. | Flexible permission system on LaminHub.
_Siloed semantics._ | External data can't be mapped on in-house data and vice versa. | Provide intuitive curation and ingestion APIs, operate on open-source data models that can be adopted by any organization, LaminHub offers overlapping schema.

## Scientific effectiveness

<!-- prettier-ignore -->
 Problem | Description | Solution
--- | --- | ---
_Support learning._ | There is no support for the learning-from-data cycle. | LaminDB supports data models across the full lab cycle, including measured → relevant → derived features. It supports managing knowledge through it's rich semantic models that interface high-dimensional data to learn from.
_Optimal decision making._ | There is no framework for optimizing decision making in the R&D organization. | LaminDB tracks all information flow including all decision making through humans and statistical models. The resulting graph models R&D operations as a whole and can in principle be optimized with sufficient data.

## Day-to-day R&D operations

<!-- prettier-ignore -->
 Problem | Description | Solution
--- | --- | ---
_Development data._ | Development data can't be ingested as data models are too rigid. | LaminDB allows ingesting development data and decorates them with corresponding QC flags to distinguish them from production data.
_Corrupted data._ | Data is corrupted | LaminDB's full provenance allows to trace back corruption to its origin and write a simple fix, typically, in form of an ingestion constraint

## Building an R&D data platform

<!-- prettier-ignore -->
 Problem | Description | Solution
--- | --- | ---
_Aligning data models._ | Data models are hard to align across interdisciplinary stakeholders. | LaminDB's templates should cover 90% of cases, the remaining 10% can be get configured.
_Lock-in._ | Commercial solutions lock organizations into specific cloud infrastructure | LaminDB works based on an open source and a multi-cloud stack with zero lock-in danger.
_Migrations are a pain._ | Migrating data models in a fast-paced R&D environment is cumbersome. | LaminDB supports schema-module-based migrations.
