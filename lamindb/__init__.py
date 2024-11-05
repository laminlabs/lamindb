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
__version__ = "0.76.15"

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
    # note regarding the ._xxx imports here:
    #   because each of these imports is dynamically changing the class interface
    #   of the django models, all of these imports actually become order dependent
    #   For that reason in this PR we have to provide Artifact in the global
    #   namespace here first, to not break consecutive imports...
    #   This would all go away when this the monkey patching would be removed from
    #   all classes.
    from . import _artifact
    Artifact = _artifact.Artifact
    from . import (
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
