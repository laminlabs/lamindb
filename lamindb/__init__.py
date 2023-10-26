"""A data framework for biology.

We assume that data is stored as files or in array formats like parquet, zarr,
HDF5, TileDB or DuckDB.

LaminDB helps you manage these data with registries for metadata.

The two most important are:

.. autosummary::
   :toctree: .

   File
   Dataset

Four registries track provenance of data batches:

.. autosummary::
   :toctree: .

   Transform
   Run
   User
   Storage

Four registries validate & contextualize data batches:

.. autosummary::
   :toctree: .

   ULabel
   Feature
   FeatureSet
   Modality

Functional tools:

.. autosummary::
   :toctree: .

   track
   view
   save

Static classes & modules:

.. autosummary::
   :toctree: .

   settings
   setup
   dev

"""

__version__ = "0.58.2"  # denote a release candidate for 0.1.0 with 0.1rc1

import os as _os

import lamindb_setup as _lamindb_setup

# prints warning of python versions
from lamin_utils import py_version_warning as _py_version_warning
from lamindb_setup import _check_instance_setup
from lamindb_setup._check_instance_setup import _INSTANCE_NOT_SETUP_WARNING
from lamindb_setup._init_instance import reload_schema_modules as _reload_schema_modules

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
        Modality,
        Run,
        Storage,
        Transform,
        ULabel,
        User,
    )

    from . import _dataset  # noqa
    from . import _feature  # noqa
    from . import _feature_set  # noqa
    from . import _file  # noqa
    from . import _parents  # noqa
    from . import _registry  # noqa
    from . import _run  # noqa
    from . import _storage  # noqa
    from . import _transform  # noqa
    from . import _ulabel  # noqa
    from . import _validate  # noqa
    from . import dev  # noqa
    from ._delete import delete  # noqa
    from ._save import save  # noqa
    from ._view import view  # noqa
    from .dev import _priors  # noqa
    from .dev._run_context import run_context  # noqa
    from .dev._settings import settings

    # schema modules
    _reload_schema_modules(_lamindb_setup.settings.instance)

    track = run_context._track  # noqa
    settings.__doc__ = """Global :class:`~lamindb.dev.Settings`."""
