"""A data framework for biology.

Core registries.

.. autosummary::
   :toctree: .

   Artifact
   Collection
   Transform
   Run
   User
   Storage
   ULabel
   Feature
   FeatureSet
   Param

Key functionality.

.. autosummary::
   :toctree: .

   track
   finish
   connect
   Curator
   view
   save

Modules and settings.

.. autosummary::
   :toctree: .

   integrations
   context
   settings
   setup
   UPath
   core

"""

# denote a release candidate for 0.1.0 with 0.1rc1, 0.1a1, 0.1b1, etc.
__version__ = "0.76.10"

import os as _os

import lamindb_setup as _lamindb_setup
from lamindb_setup._check_setup import InstanceNotSetupError as _InstanceNotSetupError
from lamindb_setup._check_setup import _check_instance_setup
from lamindb_setup._connect_instance import connect
from lamindb_setup.core.upath import UPath

from . import setup


def __getattr__(name):
    raise _InstanceNotSetupError()


if _check_instance_setup(from_module="lnschema_core"):
    del _InstanceNotSetupError
    del __getattr__  # delete so that imports work out
    from lnschema_core.models import (
        Artifact,
        Collection,
        Feature,
        FeatureSet,
        Param,
        Run,
        Storage,
        Transform,
        ULabel,
        User,
    )

    from . import core  # isort: split
    from . import (
        _artifact,
        _can_validate,
        _collection,
        _curate,
        _feature,
        _feature_set,
        _is_versioned,
        _parents,
        _record,
        _run,
        _storage,
        _transform,
        _ulabel,
        integrations,
    )
    from ._curate import Curator
    from ._save import save
    from ._view import view
    from .core._context import context
    from .core._settings import settings

    track = context.track  # simple access because these are so common
    finish = context.finish  # simple access because these are so common
    Curate = Curator  # backward compat
    settings.__doc__ = """Global settings (:class:`~lamindb.core.Settings`)."""
    context.__doc__ = """Global run context (:class:`~lamindb.core.Context`).

    Note that you can access:

    - `ln.context.track()` as `ln.track()`
    - `ln.context.finish()` as `ln.finish()`

    """
    from django.db.models import Q
