"""Schema.

Access default schema modules:

.. autosummary::
   :toctree: .

   core
   bionty
   wetlab
   bfx

Helper tools:

.. autosummary::
   :toctree: .

   view
   list_entities

"""
import lnbfx.schema as bfx

bfx.__doc__ = f"""Bioinformatics workflows.

See `lnbfx.schema <https://lamin.ai/docs/lnbfx/api>`__.
"""
import lnschema_bionty as bionty

bionty.__doc__ = f"""Biological entities.

See `lnschema-bionty <https://lamin.ai/docs/lnschema-bionty/api>`__.
"""
import lnschema_core as core

core.__doc__ = f"""Data provenance & flow.

See `lnschema-core <https://lamin.ai/docs/lnschema-core/api>`__.
"""
import lnschema_wetlab as wetlab

wetlab.__doc__ = f"""Generic wetlab.

See `lnschema-wetlab <https://lamin.ai/docs/lnschema-wetlab/api>`__.
"""

from ._core import list_entities, view

try:
    import lnschema_retro as retro
except ModuleNotFoundError:
    pass

try:
    import maren.schema as swarm
except ModuleNotFoundError:
    pass

try:
    import lnschema_harmonic_docking as docking
except ModuleNotFoundError:
    pass
