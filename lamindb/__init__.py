"""A data framework for biology.

Registries:

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

Key functionality:

.. autosummary::
   :toctree: .

   connect
   track
   finish
   Curate
   view
   save

Modules & settings:

.. autosummary::
   :toctree: .

   integrations
   settings
   setup
   UPath
   core

"""

# denote a release candidate for 0.1.0 with 0.1rc1, 0.1a1, 0.1b1, etc.
__version__ = "0.76.0"

import os as _os

import lamindb_setup as _lamindb_setup
from lamindb_setup._check_setup import InstanceNotSetupError as _InstanceNotSetupError
from lamindb_setup._check_setup import _check_instance_setup
from lamindb_setup._connect_instance import connect
from lamindb_setup.core.upath import UPath

from . import setup


def __getattr__(name):
    raise _InstanceNotSetupError()


if _check_instance_setup(from_lamindb=True):
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
    from ._curate import Curate
    from ._finish import finish
    from ._save import save
    from ._view import view
    from .core._context import context
    from .core._settings import settings

    # schema modules
    if not _os.environ.get("LAMINDB_MULTI_INSTANCE") == "true":
        from lamindb_setup._init_instance import (
            reload_schema_modules as _reload_schema_modules,
        )

        _reload_schema_modules(_lamindb_setup.settings.instance)

    track = context._track  # backward compat
    settings.__doc__ = """Global :class:`~lamindb.core.Settings`."""
    context.__doc__ = """Global :class:`~lamindb.core.Context`."""
    from django.db.models import Q
