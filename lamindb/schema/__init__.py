"""Inspect the data models for each entity.

.. currentmodule:: lamindb

For a schematic of relations among models, see :class:`track.schema`.

.. currentmodule:: lamindb.schema

See the API reference of the schema modules:

- `schema-core <https://lamin.ai/docs/lndb-schema-core/api>`__
- `schema-bionty <https://lamin.ai/docs/lndb-schema-bionty/api>`__
- `schema-wetlab <https://lamin.ai/docs/lndb-schema-wetlab/api>`__
- `schema-bfx <https://lamin.ai/docs/lndb_bfx_pipeline/api>`__

"""
import lndb_bfx_pipeline.schema as bfx
import lndb_schema_bionty as bionty
import lndb_schema_core as core
import lndb_schema_wetlab as wetlab
