"""Schema tools & overview.

Guide: :doc:`/schemas`

You can access mounted schema modules with domain-specific entities via
available via `ln.schema.<module>.<entity>`.

However, we recommend to import schema modules, e.g., like `import
lnschema_bionty as bt`.

.. autosummary::
   :toctree: .

   graph
   view

"""
import importlib as _importlib

from lamindb_setup import settings as _settings
from lamindb_setup._init_instance import reload_schema_modules as _reload_schema_modules

from .. import _INSTANCE_SETUP

if _INSTANCE_SETUP:
    _reload_schema_modules(_settings.instance)

from ._core import graph, view
