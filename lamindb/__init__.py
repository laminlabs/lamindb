"""LaminDB: Manage R&D data & analyses.

Import the package::

   import lamindb as ln

`File` tracks data artifacts in form of files, on-disk (`zarr`, etc.) and
in-memory data objects (`DataFrame`, `AnnData`, etc.) and allows to link them
against entities of core schema & custom schemas.

The core schema entities are central to lamindb's API:

.. autosummary::
   :toctree: .

   File
   Dataset
   Transform
   Run
   Feature
   Label
   User
   Storage

More control over feature management:

.. autosummary::
   :toctree: .

   FeatureSet

Functional tools:

.. autosummary::
   :toctree: .

   track
   view
   save
   delete

Static classes & modules:

.. autosummary::
   :toctree: .

   settings
   context
   types
   setup
   schema
   dev

"""

__version__ = "0.48a3"  # denote a release candidate for 0.1.0 with 0.1rc1

import os as _os

import lamindb_setup as _lamindb_setup

# prints warning of python versions
from lamin_utils import py_version_warning as _py_version_warning
from lamindb_setup import _check_instance_setup
from lamindb_setup._check_instance_setup import _INSTANCE_NOT_SETUP_WARNING

_py_version_warning("3.8", "3.11")

_TESTING = _lamindb_setup._TESTING
_INSTANCE_SETUP = _check_instance_setup(from_lamindb=True)
# allow the user to call setup
from . import setup  # noqa


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
    from lnschema_core import (  # noqa
        Dataset,
        Feature,
        FeatureSet,
        File,
        Label,
        Run,
        Storage,
        Transform,
        User,
    )

    from . import dev  # noqa
    from . import schema  # noqa
    from . import types  # noqa
    from ._context import context  # noqa

    track = context._track  # noqa
    from lamin_utils import logger as _logger

    from . import _dataset  # noqa
    from . import _feature  # noqa
    from . import _feature_set  # noqa
    from . import _file  # noqa
    from . import _label  # noqa
    from . import _orm  # noqa
    from . import _transform  # noqa
    from ._delete import delete  # noqa
    from ._save import save  # noqa
    from ._select import select  # noqa
    from ._view import view  # noqa
    from .dev._settings import settings

    add = save  # backward compat

    settings.__doc__ = """Global :class:`~lamindb.dev.Settings`."""
