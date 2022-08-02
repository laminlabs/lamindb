# Ingestion

**Problem:** We have some data from somewhere and we want to store it. We call this process _ingestion_.

**Context:** External data that is subject to external conventions (schema & formats) and internal data (generation & pipelining constrained to internal schema & formats) may need ingestion processes with highly differing levels of _curation_: cleaning, standardization, normalization, and annotation.

**Examples:** In this tutorial, we'll choose examples on the whole spectrum!

Let us start to solve this problem with a local database and storage.

A few defintions & rationales:

- Data curation is the process that makes data ingestion succeed.
- If you don't curate sufficiently (and you can configure what the means), you'll not be able to ingest.
