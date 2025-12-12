"""A data framework for biology.

If you just want to read data from a LaminDB instance, use QueryDB, no setup required::

   import lamindb as ln

   db = ln.QueryDB("laminlabs/cellxgene")

Here is the reference.

.. autosummary::
   :toctree: .

   QueryDB

If you want to write data to a LaminDB instance, you need collaborator access, run:

```shell
lamin login
lamin connect account/name
```

You can create an instance at [lamin.ai](https://lamin.ai) similar to how you create repo on GitHub.
If you just want to work with a local SQLite instance, run:

```shell
lamin init --storage ./quickstart-data --modules bionty
```


Lineage
=======

Track inputs, outputs & environment of a notebook or script run.

.. autosummary::
   :toctree: .

   track
   finish

Decorate a function with `@tracked()` to track inputs, outputs & environment of function executions.

.. autosummary::
   :toctree: .

   tracked

Artifacts & storage locations
=============================

Files, folders & arrays and their storage locations.

.. autosummary::
   :toctree: .

   Artifact
   Storage

Transforms & runs
=================

Data transformations and their executions.

.. autosummary::
   :toctree: .

   Transform
   Run

Records, labels, features & schemas
===================================

Create labels and manage flexible records, e.g., for samples or donors.

.. autosummary::
   :toctree: .

   Record
   ULabel

Define features & schemas to validate artifacts & records.

.. autosummary::
   :toctree: .

   Feature
   Schema

Project management
==================

.. autosummary::
   :toctree: .

   User
   Collection
   Project
   Space
   Branch
   Reference

Basic utilities
===============

Connecting, viewing database content, accessing settings & run context.

.. autosummary::
   :toctree: .

   connect
   view
   save
   UPath
   settings
   context

Curators and integrations
=========================

.. autosummary::
   :toctree: .

   curators
   integrations

Examples, errors & setup
========================

.. autosummary::
   :toctree: .

   examples
   errors
   setup

Developer API
=============

.. autosummary::
   :toctree: .

   base
   core
   models

"""

# ruff: noqa: I001
# denote a release candidate for 0.1.0 with 0.1rc1, 0.1a1, 0.1b1, etc.
__version__ = "1.17a1"

import warnings as _warnings

# through SpatialData
_warnings.filterwarnings(
    "ignore", message="The legacy Dask DataFrame implementation is deprecated"
)

from lamindb_setup._check_setup import _check_instance_setup
from lamindb_setup._connect_instance import connect
from lamindb_setup.core.upath import UPath

from . import base, errors, setup

_check_instance_setup(from_module="lamindb")

from ._tracked import tracked
from ._view import view
from .core._context import context
from .core._settings import settings
from .curators._legacy import CatManager as Curator
from .models import (
    Artifact,
    Collection,
    Feature,
    FeatureSet,  # backward compat
    Project,
    Reference,
    Run,
    Schema,
    Storage,
    Transform,
    ULabel,
    User,
    Space,
    Branch,
    Record,
    QueryDB,
)
from .models.save import save
from . import core
from . import integrations
from . import curators
from . import examples

track = context._track
finish = context._finish
settings.__doc__ = """Global live settings (:class:`~lamindb.core.Settings`)."""
context.__doc__ = """Global run context (:class:`~lamindb.core.Context`)."""

from django.db.models import Q

Param = Feature  # backward compat

__all__ = [
    # data lineage
    "track",
    "finish",
    "tracked",
    # registries
    "Artifact",
    "Storage",
    "Transform",
    "Run",
    "Feature",
    "ULabel",
    "Schema",
    "Record",
    "User",
    "Collection",
    "Project",
    "Space",
    "Branch",
    "Reference",
    # other
    "connect",
    "view",
    "save",
    "UPath",
    "settings",
    "context",
    "QueryDB",
    # curators and integrations
    "curators",
    "integrations",
    # examples, errors, setup
    "examples",
    "errors",
    "setup",
    # low-level functionality
    "base",
    "core",
    "models",
]
