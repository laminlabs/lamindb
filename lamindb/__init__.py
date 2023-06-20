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
   Transform
   Run
   User
   Storage
   FeatureSet
   Folder
   Project

Functional tools:

.. autosummary::
   :toctree: .

   track
   parse
   view
   select
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

Developer API:

.. autosummary::
   :toctree: .

   dev

"""

__version__ = "0.44.0"  # denote a release candidate for 0.1.0 with 0.1rc1

import lamindb_setup as _lamindb_setup

# prints warning of python versions
from lamin_logger import py_version_warning as _py_version_warning
from lamindb_setup import _check_instance_setup

_py_version_warning("3.8", "3.10")

_INSTANCE_SETUP = _check_instance_setup(from_lamindb=True)
# allow the user to call setup
from . import setup  # noqa

# only import all other functionality if setup was successful
if _INSTANCE_SETUP:
    from lnschema_core import (  # noqa
        FeatureSet,
        File,
        Folder,
        Project,
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
    from lamin_logger import logger as _logger

    # this needs to follow on the import right now
    _logger.success(
        f"Loaded instance: {_lamindb_setup.settings.instance.identifier} (lamindb"
        f" {__version__})"
    )

    from . import _baseorm_methods  # noqa
    from . import _featureset_methods  # noqa
    from . import _file_methods  # noqa
    from . import _folder_methods  # noqa
    from . import _transform_methods  # noqa
    from ._delete import delete  # noqa
    from ._parse import parse  # noqa
    from ._save import save  # noqa
    from ._select import select  # noqa
    from ._view import view  # noqa
    from .dev._settings import settings

    add = save  # backward compat

    settings.__doc__ = """Global :class:`~lamindb.dev.Settings`."""
