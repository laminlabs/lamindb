"""Inspect the data models for each entity.

.. autosummary::
   :toctree: .

   draw
   list_entities

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

from ._core import draw, list_entities

try:
    import lnschema_retro as retro
except ModuleNotFoundError:
    pass

try:
    import maren.schema as swarm
except ModuleNotFoundError:
    pass
