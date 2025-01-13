"""A data framework for biology.

Tracking notebooks & scripts.

.. autosummary::
   :toctree: .

   track
   finish

Registries.

.. autosummary::
   :toctree: .

   Artifact
   Transform
   ULabel
   Run
   User
   Storage
   Feature
   Schema
   Param
   Collection
   Project
   Reference
   Person

Key functionality.

.. autosummary::
   :toctree: .

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
   base
   core

"""

# denote a release candidate for 0.1.0 with 0.1rc1, 0.1a1, 0.1b1, etc.
__version__ = "1.0a4"

from lamindb_setup._check_setup import InstanceNotSetupError as _InstanceNotSetupError
from lamindb_setup._check_setup import _check_instance_setup
from lamindb_setup._connect_instance import connect
from lamindb_setup.core.upath import UPath

from . import base, setup


def __getattr__(name):
    raise _InstanceNotSetupError()


if _check_instance_setup(from_module="lamindb"):
    del __getattr__  # so that imports work out
    from . import core  # isort: split
    from . import (
        _artifact,
        _can_curate,
        _collection,
        _feature,
        _is_versioned,
        _parents,
        _record,
        _run,
        _schema,
        _storage,
        _transform,
        _ulabel,
        integrations,
    )
    from ._save import save
    from ._view import view
    from .core._context import context
    from .core._settings import settings
    from .curators import Curator
    from .models import (
        Artifact,
        Collection,
        Feature,
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
    )

    track = context.track  # simple access
    finish = context.finish  # simple access
    settings.__doc__ = """Global settings (:class:`~lamindb.core.Settings`)."""
    context.__doc__ = """Global run context (:class:`~lamindb.core.Context`).

    Note that you can access:

    - `ln.context.track()` as `ln.track()`
    - `ln.context.finish()` as `ln.finish()`

    """
    from django.db.models import Q
