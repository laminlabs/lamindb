"""Schema.

Import the submodule::

   import lamindb.schema as lns

- The entities in the `core schema <https://lamin.ai/docs/lnschema-core>`__ are
  the basis for tracking any data and are available from `lns.<entity>`.
- Additional mounted schema modules provide domain-specific entities and are
  available via `lns.<module>.<entity>`.

Data objects & transformations:

.. autosummary::
   :toctree: .

   DSet
   DObject
   DTransform
   DTransformIn

Users, projects, storage locations, and usage statistics:

.. autosummary::
   :toctree: .

   User
   Project
   Storage
   Usage

Default data transformations:

.. autosummary::
   :toctree: .

   Jupynb
   Pipeline
   PipelineRun

Any LaminDB schema module that has been mounted to an instance can be accessed like the bionty, wetlab, bfx modules below:

.. autosummary::
   :toctree: .

   bionty
   wetlab
   bfx

Helper tools:

.. autosummary::
   :toctree: .

   view
   list_tables
   dev

"""
import lnbfx.schema as bfx
from lnschema_core import (
    DObject,
    DSet,
    DSetDObject,
    DTransform,
    DTransformIn,
    Jupynb,
    Pipeline,
    PipelineRun,
    Project,
    ProjectDSet,
    Storage,
    Usage,
    User,
    dev,
)

bfx.__doc__ = f"""Bioinformatics workflows.

See `lnbfx.schema <https://lamin.ai/docs/lnbfx/api>`__.
"""
import lnschema_bionty as bionty

bionty.__doc__ = f"""Biological entities.

See `lnschema-bionty <https://lamin.ai/docs/lnschema-bionty/api>`__.
"""

import lnschema_wetlab as wetlab

wetlab.__doc__ = f"""Generic wetlab.

See `lnschema-wetlab <https://lamin.ai/docs/lnschema-wetlab/api>`__.
"""

from ._core import list_tables, view

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


list_entities = list_tables  # backward compat
