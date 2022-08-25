"""Inspect the data models for each entity.

.. currentmodule:: lamindb

For a schematic of relations among models, see :class:`track.schema`.

.. currentmodule:: lamindb.schema

See the API reference of the schema modules:

- `schema-core <https://lamin.ai/docs/lnschema-core/api>`__
- `schema-bionty <https://lamin.ai/docs/lnschema-bionty/api>`__
- `schema-wetlab <https://lamin.ai/docs/lnschema-wetlab/api>`__
- `schema-bfx <https://lamin.ai/docs/lnbfx/api>`__

"""
import lnbfx.schema as bfx
import lnschema_bionty as bionty
import lnschema_core as core
import lnschema_wetlab as wetlab

try:
    import lnschema_retro as retro
except ModuleNotFoundError:
    pass
