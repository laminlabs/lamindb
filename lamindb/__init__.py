"""A data framework for biology.

Installation::

   pip install lamindb

If you just want to *read* data from a LaminDB instance, use :class:`~lamindb.DB`::

   import lamindb as ln

   db = ln.DB("laminlabs/cellxgene")

To *write* data, connect to a writable instance::

   lamin login
   lamin connect account/name

You can create an instance at `lamin.ai <https://lamin.ai>`__ and invite collaborators.
If you prefer to work with a local database (no login required), run::

    lamin init --storage ./quickstart-data --modules bionty

LaminDB will then auto-connect upon import and you can then create & save objects like this::

   import lamindb as ln
   # → connected lamindb: account/instance

   ln.Artifact("./my_dataset.parquet", key="datasets/my_dataset.parquet").save()

Lineage
=======

Track inputs, outputs, parameters, and environments of notebooks, scripts, and functions.

.. autosummary::
   :toctree: .

   track
   finish
   flow
   step

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

   DB
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
__version__ = "2.1.0"

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

from .core._functions import flow, step, tracked
from ._view import view
from .core._context import context
from .core._settings import settings
from .models import (
    Artifact,
    Collection,
    Feature,
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
    DB,
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
    "step",
    "flow",
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
    "DB",
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
