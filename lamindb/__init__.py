"""A data framework for biology.

Tracking notebooks, scripts & functions.

.. autosummary::
   :toctree: .

   track
   finish
   tracked

Registries.

.. autosummary::
   :toctree: .

   Artifact
   Transform
   Run
   ULabel
   User
   Storage
   Feature
   Schema
   Param
   Collection
   Project
   Space
   Reference
   Person

Curators & integrations.

.. autosummary::

   curators
   integrations

Key functionality.

.. autosummary::
   :toctree: .

   connect
   view
   save
   UPath
   settings

Low-level functionality.

.. autosummary::
   :toctree: .

   context
   errors
   setup
   base
   models
   core

Backward compatibility.

.. autosummary::
   :toctree: .

   FeatureSet
   Curator

"""

# ruff: noqa: I001
# denote a release candidate for 0.1.0 with 0.1rc1, 0.1a1, 0.1b1, etc.
__version__ = "1.2a2"

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
    from .curators import CatManager as Curator
    from .models import (
        Artifact,
        Collection,
        Feature,
        FeatureSet,  # backward compat
        Param,
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
    )
    from .models.save import save
    from . import core
    from . import integrations

    track = context.track  # simple access
    finish = context.finish  # simple access
    settings.__doc__ = """Global settings (:class:`~lamindb.core.Settings`)."""
    context.__doc__ = """Global run context (:class:`~lamindb.core.Context`).

    Note that you can access:

    - `ln.context.track()` as `ln.track()`
    - `ln.context.finish()` as `ln.finish()`

    """
    from django.db.models import Q
