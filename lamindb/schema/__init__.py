"""Access the schema.

- The core entities are the basis for tracking any data and are available from
  `lns.<entity>`.
- Additional mounted schema modules provide domain-specific entities and are
  available via `lns.<module>.<entity>`.

Import this submodule as::

   import lamindb.schema as lns

Core entities
=============

Data objects ({class}`lamindb.DObject`) are transformed by runs:

.. autosummary::
   :toctree: .

   Run

Runs transform data using code:

.. autosummary::
   :toctree: .

   Pipeline
   Notebook

Grouping data objects as sets or by features:

.. autosummary::
   :toctree: .

   DSet
   Features

Users, projects, storage locations, and usage:

.. autosummary::
   :toctree: .

   User
   Project
   Storage
   Usage

See the source code `here <https://lamin.ai/docs/lnschema-core>`__.

Exemplary extensions
====================

Any LaminDB schema module that has been mounted to an instance can be accessed like the bionty, wetlab, bfx modules below:

- `bionty <https://lamin.ai/docs/lnschema-bionty/api>`__: Knowledge-managed biological entities.
- `wetlab <https://lamin.ai/docs/lnschema-wetlab/api>`__: Wetlab operations.


Helper tools
============

.. autosummary::
   :toctree: .

   view
   list_tables
   dev

"""
import importlib as _importlib

from lndb_setup import settings as _settings
from lndb_setup._setup_schema import get_schema_module_name as _get_schema_module_name
from lnschema_core import (
    DObject,
    DSet,
    Features,
    Notebook,
    Pipeline,
    Project,
    Run,
    Storage,
    Usage,
    User,
    dev,
)
from packaging import version as _v

_check_v = {
    "bionty": "0.6.8",
    "wetlab": "0.13.3",
    "bfx": "0.7.1",
}

for name in _settings.instance.schema:
    _module = _importlib.import_module(_get_schema_module_name(name))
    if name in _check_v:
        if _v.parse(_module.__version__) != _v.parse(_check_v[name]):
            raise RuntimeError(f"lamindb needs lnschema_{name}=={_check_v[name]}")
    globals()[_module._name] = _module

from ._core import list_tables, view

list_entities = list_tables  # backward compat
