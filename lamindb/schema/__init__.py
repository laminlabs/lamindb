"""Schema tools & overview.

Guide: :doc:`/guide/schema`

You can access mounted schema modules with domain-specific entities via
available via `ln.schema.<module>.<entity>`.

However, we recommend to import schema modules, e.g., like `import
lnschema_bionty as bt`.

.. autosummary::
   :toctree: .

   view
   dev

"""
import importlib as _importlib

from lndb import settings as _settings
from lndb.dev._setup_schema import (
    check_schema_version_and_import as _check_schema_version_and_import,
)
from lnschema_core import Features, Project, Run, Storage
from lnschema_core import Transform as _Transform
from lnschema_core import User, dev

from .. import _instance_setup


def _import_schema():
    for name in _settings.instance.schema:
        module = _check_schema_version_and_import(name)
        globals()[module._name] = module


if _instance_setup:
    _import_schema()

from ._core import list_tables, view

list_entities = list_tables  # backward compat

Pipeline = _Transform
Notebook = _Transform
