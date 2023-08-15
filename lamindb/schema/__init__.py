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
from lamindb_setup.dev._setup_schema import (
    check_schema_version_and_import as _check_schema_version_and_import,
)

from .. import _INSTANCE_SETUP


def _import_schema():
    for name in _settings.instance.schema:
        module = _check_schema_version_and_import(name)
        globals()[module._name] = module


if _INSTANCE_SETUP:
    _import_schema()

from ._core import graph, view
