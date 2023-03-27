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
from lndb.dev._setup_schema import get_schema_module_name as _get_schema_module_name
from lnschema_core import Features, Project, Run, Storage
from lnschema_core import Transform as _Transform
from lnschema_core import User, dev
from packaging import version as _v

_check_v = {
    "bionty": "0.13.1",
    "wetlab": "0.15rc2",
}

from .. import _instance_setup


def _import_schema():
    for name in _settings.instance.schema:
        _module = _importlib.import_module(_get_schema_module_name(name))
        if name in _check_v:
            if _v.parse(_module.__version__) != _v.parse(_check_v[name]):
                raise RuntimeError(f"lamindb needs lnschema_{name}=={_check_v[name]}")
        globals()[_module._name] = _module


if _instance_setup:
    _import_schema()

from ._core import list_tables, view

list_entities = list_tables  # backward compat

Pipeline = _Transform
Notebook = _Transform
