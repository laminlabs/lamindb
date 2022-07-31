# Common problems of data-heavy R&D organizations

We'll add the summary here once we completed the detailed table.

## Detailed table

<!-- prettier-ignore -->
| | LaminDB's solution
--- | ---
**Fundamental scientific problems** |
Key results aren't supported by solid evidence: analytics results (p-values, effectiveness, etc.) cannot be linked to supporting data as too many processing steps and stakeholders are involved | LaminDB provides full provenance
Key data can't be accessed (queried) as it's buried in blob storage | Index observations and variables in blob storage and link them against a query database
Key data can't be accessed as there is no guiding structure | Structure data both through the scientific  hypothesis space and through provenance (data flow and sequential R&D operations)
Data can't be accessed at scale as no viable programmatic interfaces exist | LaminDB creates an API basis that can power programmatic access and graphical wetlab interfaces
Molecular (high-dimensional) data can't be integrated with phenotypic (low-dimensional) data | Index molecular data with the same entities as phenotypic data, provide interfaces with databases for low-dimensional data
**Problems in collaborative science** |
Data can't be easily shared across organizations in an interoperable and reusable way | Flexible permission system on LaminHub, data is encoded with widely overlapping schema
**Problems in day-to-day science** |
Data can't be queried by relevant entities | Provide rich default data models with a configurable template system, in particular, provide a data model for biological entities
External data can't be mapped on in-house data | Provide intuitive curation and ingestion APIs, operate on open-source data models that can be adopted by any organization
There is no support for the typical learning-from-data cycle | LaminDB supports data models across the full lab cycle, including measured → relevant → derived features, and supports expanding knowledge management through it's rich semantic models along side managing high-dimensional data to learn from
**Problems in day-to-day R&D operations** |
Development data can't be ingested as data models are too rigid | LaminDB allows ingesting "draft" data and decorates them with corresponding QC flags to distinguish them from "production data"
Data is corrupted | LaminDB's full provenance allows to trace back corruption to its origin and write a simple fix, typically, in form of an ingestion constraint
**Problems in building data platforms** |
Data models are hard to align as data models require highly cross-functional input | LaminDB's templates should cover 90% of cases, the remaining 10% can be get configured
Commercial solutions lock organizations into specific cloud infrastructure | LaminDB works based on an open source and a multi-cloud stack with zero lock-in danger
