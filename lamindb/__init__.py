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
   dev

"""

__version__ = "0.67.0"  # denote a release candidate for 0.1.0 with 0.1rc1

import os as _os

import lamindb_setup as _lamindb_setup

# prints warning of python versions
from lamin_utils import py_version_warning as _py_version_warning
from lamindb_setup import _check_instance_setup
from lamindb_setup._check_instance_setup import _INSTANCE_NOT_SETUP_WARNING
from lamindb_setup._init_instance import reload_schema_modules as _reload_schema_modules
from lamindb_setup.dev.upath import UPath

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
    from lnschema_core import (
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

    File = Artifact  # backward compat
    from . import _artifact  # noqa
    from . import _collection
    from . import _feature
    from . import _feature_set
    from . import _parents
    from . import _registry
    from . import _run
    from . import _storage
    from . import _transform
    from . import _ulabel
    from . import _validate
    from . import dev
    from ._delete import delete
    from ._save import save
    from ._view import view
    from .dev._run_context import run_context
    from .dev._settings import settings

    # schema modules
    _reload_schema_modules(_lamindb_setup.settings.instance)

    track = run_context._track
    settings.__doc__ = """Global :class:`~lamindb.dev.Settings`."""
