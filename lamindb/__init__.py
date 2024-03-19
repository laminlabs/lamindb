"""A data framework for biology.

Registries track artifacts & metadata across provenance, features & labels.

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

Functions:

.. autosummary::
   :toctree: .

   connect
   track
   finish
   view
   save

Modules & settings:

.. autosummary::
   :toctree: .

   validation
   settings
   setup
   UPath
   core

"""

__version__ = "0.69.1"  # denote a release candidate for 0.1.0 with 0.1rc1

import os as _os

import lamindb_setup as _lamindb_setup
from lamindb_setup._check_setup import InstanceNotSetupError, _check_instance_setup
from lamindb_setup._connect_instance import connect
from lamindb_setup.core.upath import UPath

from . import setup


def __getattr__(name):
    raise InstanceNotSetupError()


if _check_instance_setup(from_lamindb=True):
    del InstanceNotSetupError
    del __getattr__  # delete so that imports work out
    from lnschema_core.models import (
        Artifact,
        Collection,
        Feature,
        FeatureSet,
        Run,
        Storage,
        Transform,
        ULabel,
        User,
    )

    from . import (
        _artifact,
        _collection,
        _feature,
        _feature_set,
        _is_versioned,
        _parents,
        _registry,
        _run,
        _storage,
        _transform,
        _ulabel,
        _validate,
        core,
    )

    dev = core  # backward compat
    from ._finish import finish
    from ._save import save
    from ._view import view
    from .core._run_context import run_context as _run_context
    from .core._settings import settings
    from .core._transform_settings import transform  # backward compat

    # schema modules
    if not _os.environ.get("LAMINDB_MULTI_INSTANCE") == "true":
        from lamindb_setup._init_instance import (
            reload_schema_modules as _reload_schema_modules,
        )

        _reload_schema_modules(_lamindb_setup.settings.instance)

    track = _run_context._track
    settings.__doc__ = """Global :class:`~lamindb.core.Settings`."""
