"""A data framework for biology.

Data lineage
============

Track inputs, outputs & environment of a notebook or script run.

.. autosummary::
   :toctree: .

   track
   finish

Decorate a function with `@tracked()` to track inputs, outputs & environment of function executions.

.. autosummary::
   :toctree: .

   tracked

Registries
==========

Manage artifacts and transforms.

.. autosummary::
   :toctree: .

   Artifact
   Storage
   Transform
   Run

Validate and annotate artifacts.

.. autosummary::
   :toctree: .

   Feature
   ULabel
   Schema

Manage flexible records to track, e.g., samples or donors.

.. autosummary::
   :toctree: .

   Record
   Sheet

Manage projects.

.. autosummary::
   :toctree: .

   User
   Collection
   Project
   Space
   Branch
   Reference
   Person

Other
=====

Functions and classes.

.. autosummary::
   :toctree: .

   connect
   view
   save
   UPath
   settings
   context

Curators and integrations.

.. autosummary::
   :toctree: .

   curators
   integrations

Low-level functionality.

.. autosummary::
   :toctree: .

   examples
   errors
   setup
   base
   core
   models

Backwards compatibility.

.. autosummary::
   :toctree: .

   Param
   FeatureSet
   Curator

"""

# ruff: noqa: I001
# denote a release candidate for 0.1.0 with 0.1rc1, 0.1a1, 0.1b1, etc.
__version__ = "1.6.1"

import warnings

# through SpatialData
warnings.filterwarnings(
    "ignore", message="The legacy Dask DataFrame implementation is deprecated"
)

from lamindb_setup._check_setup import InstanceNotSetupError as _InstanceNotSetupError
from lamindb_setup._check_setup import _check_instance_setup
from lamindb_setup._connect_instance import connect
from lamindb_setup.core.upath import UPath

from . import base, errors, setup


def __getattr__(name):
    raise _InstanceNotSetupError()


if _check_instance_setup(from_module="lamindb"):
    del __getattr__  # so that imports work out
    from . import base
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
        Person,
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
        Sheet,
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
