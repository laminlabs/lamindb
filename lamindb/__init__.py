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

from __future__ import annotations

import importlib
import warnings
from typing import TYPE_CHECKING

from lamindb_setup import connect
from lamindb_setup._check_setup import _check_instance_setup
from lamindb_setup.core.upath import UPath

# denote a release candidate for 0.1.0 with 0.1rc1, 0.1a1, 0.1b1, etc.
__version__ = "1.6.1"

_loaded_submodules_pre_connect = [
    "base",
    "errors",
    "setup",
]
_loaded_submodules_post_connect = [
    "core",
    "integrations",
    "curators",
    "examples",
]

__all__ = [
    # available without instance connection
    "connect",
    "UPath",
    # models
    "Artifact",
    "Collection",
    "Feature",
    "Person",
    "Project",
    "Reference",
    "Run",
    "Schema",
    "Storage",
    "Transform",
    "ULabel",
    "User",
    "Space",
    "Branch",
    "Record",
    "Sheet",
    # models django.db
    "Q",
    # functions
    "tracked",
    "view",
    "save",
    "track",
    "finish",
    # objects
    "context",
    "settings",
    # backwards compat
    "FeatureSet",
    "Param",
    "Curator",
    # extras
    *_loaded_submodules_pre_connect,
    *_loaded_submodules_post_connect,
]


def __dir__():
    return __all__


# through SpatialData
warnings.filterwarnings(
    "ignore", message="The legacy Dask DataFrame implementation is deprecated"
)


def __getattr__(name):
    if name in _loaded_submodules_pre_connect:
        globals()[name] = mod = importlib.import_module(f"lamindb.{name}")
        return mod

    elif name not in __all__:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    else:
        from lamindb_setup._check_setup import InstanceNotSetupError

        raise InstanceNotSetupError()


if _check_instance_setup(from_module="lamindb") or TYPE_CHECKING:
    import lamindb.base as base
    import lamindb.core as core
    import lamindb.curators as curators
    import lamindb.examples as examples
    import lamindb.integrations as integrations

    from ._tracked import tracked
    from ._view import view
    from .core._context import context
    from .core._settings import settings
    from .curators._legacy import CatManager as Curator
    from .models import (
        Artifact,
        Branch,
        Collection,
        Feature,
        FeatureSet,  # backward compat
        Person,
        Project,
        Record,
        Reference,
        Run,
        Schema,
        Sheet,
        Space,
        Storage,
        Transform,
        ULabel,
        User,
    )
    from .models.save import save

    track = context._track
    finish = context._finish
    settings.__doc__ = """Global live settings (:class:`~lamindb.core.Settings`)."""
    context.__doc__ = """Global run context (:class:`~lamindb.core.Context`)."""
    from django.db.models import Q

    Param = Feature  # backward compat

else:
    import lamindb._lazy_import as _lazy_import

    # classes
    for cls_name in [
        "Artifact",
        "Collection",
        "Feature",
        "Person",
        "Project",
        "Reference",
        "Run",
        "Schema",
        "Storage",
        "Transform",
        "ULabel",
        "User",
        "Space",
        "Branch",
        "Record",
        "Sheet",
        "Q",
        "FeatureSet",
        "Param",
        "Curator",
    ]:
        # classes have to be assigned dynamically to prevent mypy erroring on type assignment
        globals()[cls_name] = _lazy_import.new_metaclass_proxy(cls_name)

    # functions
    tracked = _lazy_import.new_callable_proxy("tracked")
    view = _lazy_import.new_callable_proxy("view")
    save = _lazy_import.new_callable_proxy("save")
    track = _lazy_import.new_callable_proxy("track")
    finish = _lazy_import.new_callable_proxy("finish")

    # objects
    context = _lazy_import.new_instance_proxy("context")
    settings = _lazy_import.new_instance_proxy("settings")

    del _lazy_import
