"""A data framework for biology.

LaminDB helps you manage data batches with two basic registries:

.. autosummary::
   :toctree: .

   Artifact
   Collection

Four registries track provenance of data batches:

.. autosummary::
   :toctree: .

   Transform
   Run
   User
   Storage

Three registries validate & contextualize:

.. autosummary::
   :toctree: .

   ULabel
   Feature
   FeatureSet

You can also access data directly via paths:

.. autosummary::
   :toctree: .

   UPath

Functions:

.. autosummary::
   :toctree: .

   track
   view
   save

Modules & settings:

.. autosummary::
   :toctree: .

   settings
   setup
   core

"""

__version__ = "0.68.0"  # denote a release candidate for 0.1.0 with 0.1rc1

import os as _os

import lamindb_setup as _lamindb_setup

# prints warning of python versions
from lamin_utils import py_version_warning as _py_version_warning
from lamindb_setup import _check_instance_setup
from lamindb_setup._check_instance_setup import _INSTANCE_NOT_SETUP_WARNING
from lamindb_setup._init_instance import reload_schema_modules as _reload_schema_modules
from lamindb_setup.core.upath import UPath

_py_version_warning("3.8", "3.11")

_TESTING = _lamindb_setup._TESTING
_INSTANCE_SETUP = _check_instance_setup(from_lamindb=True)

# allow the user to call setup
from . import setup


class InstanceNotSetupError(Exception):
    pass


def __getattr__(name):
    raise InstanceNotSetupError(
        f"{_INSTANCE_NOT_SETUP_WARNING}If you used the CLI to init or load an instance,"
        " please RESTART the python session (in a notebook, restart kernel)"
    )


# only import all other functionality if setup was successful
if _INSTANCE_SETUP:
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
    from ._save import save
    from ._view import view
    from .core._run_context import run_context
    from .core._settings import settings
    from .core._transform_settings import transform_settings as transform

    # schema modules
    if not _os.environ.get("LAMINDB_MULTI_INSTANCE") == "true":
        _reload_schema_modules(_lamindb_setup.settings.instance)

    track = run_context._track
    settings.__doc__ = """Global :class:`~lamindb.core.Settings`."""
