# Common problems of data-heavy R&D

## Fundamental scientific problems

<!-- prettier-ignore -->
Problem | Description | Solution
--- | --- | ---
_Stand on solid ground._ | Key analytics results cannot be linked to supporting data as too many processing steps are involved. | Provide full provenance.
_Blob storage._ | Data can't be queried as it's buried in blob storage. | Index observations and variables in blob storage and link them against a query database.
_Pile of data._ | Data can't be accessed as it's not well-structured. | Structure data both via the scientific hypothesis space and via provenance (data flow and iterative R&D operations).
_Learn from data at scale._ | Data can't be accessed at scale as no viable programmatic interfaces exist. | Rich API.
_Data integration._ | Molecular (high-dimensional) data can't be integrated with phenotypic (low-dimensional) data. | Index molecular data with the same entities as phenotypic data. Provide connectors with storage of low-dimensional data (ELN & LIMS systems).

## Collaborative science

<!-- prettier-ignore -->
Problem | Description | Solution
--- | --- | ---
_Siloed infrastructure._ | Data can't be easily shared across organizations. | Lamin Hub.
_Siloed semantics._ | External data can't be mapped on in-house data and vice versa. | Provide curation and ingestion API, operate on open-source data models that can be adopted by any organization, allow unions and joins on overlapping schema through Hub.

## Scientific effectiveness

<!-- prettier-ignore -->
Problem | Description | Solution
--- | --- | ---
_Support learning._ | There is no support for the learning-from-data cycle. | Support data models across the full lab cycle, including measured → relevant → derived features. Manage knowledge through rich semantic models that map high-dimensional data to learn from.
_Optimal decision making._ | There is no framework for optimizing decision making. | LaminDB tracks all information flow including all decision making through humans and statistical models. The resulting graph models R&D operations as a whole and can in principle be optimized with sufficient data.

## Day-to-day R&D operations

<!-- prettier-ignore -->
Problem | Description | Solution
--- | --- | ---
_Development data._ | Data associated with assay development can't be ingested as data models are too rigid. | Ingest data of any curation level and label them with corresponding QC flags.
_Corrupted data._ | Data is corrupted. | Full provenance allows to trace back corruption to its origin and write a simple fix, typically, in form of an ingestion constraint.

## Building an R&D data platform

<!-- prettier-ignore -->
Problem | Description | Solution
--- | --- | ---
_Aligning data models._ | Data models are hard to align across interdisciplinary stakeholders. | Templates cover 90% of cases, the remaining 10% can be get configured.
_Lock-in._ | Commercial solutions lock organizations into specific cloud infrastructure. | Open-source and a multi-cloud stack with zero lock-in danger.
_Migrations are a pain._ | Migrating data models in a fast-paced R&D environment is often prohibitive. | Tools for schema module migrations.
