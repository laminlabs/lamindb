"""LaminDB: Manage R&D data & analyses.

Import the package::

   import lamindb as ln

The central class of the API is `File`, a wrapper for files, on-disk (`zarr`, etc.)
and in-memory data objects (`DataFrame`, `AnnData`, etc.).

.. autosummary::
   :toctree: .

   File
   Folder

Track runs of data transformations:

.. autosummary::
   :toctree: .

   Run
   Transform

Manage feature sets:

.. autosummary::
   :toctree: .

   Featureset

Query & manipulate data:

.. autosummary::
   :toctree: .

   select
   save
   delete

Utility functions:

.. autosummary::
   :toctree: .

   parse
   track
   view

Basic entities:

.. autosummary::
   :toctree: .

   User
   Project
   Storage

Schema - entities and their relations:

.. autosummary::
   :toctree: .

   schema

Setup:

.. autosummary::
   :toctree: .

   setup

Developer API:

.. autosummary::
   :toctree: .

   context
   settings
   types
   dev
"""

__version__ = "0.42a8"  # denote a release candidate for 0.1.0 with 0.1rc1

import lamindb_setup as _lamindb_setup

# prints warning of python versions
from lamin_logger import py_version_warning as _py_version_warning
from lamindb_setup._check_instance_setup import (
    check_instance_setup as _check_instance_setup,
)

_py_version_warning("3.8", "3.10")

_INSTANCE_SETUP = _check_instance_setup(from_lamindb=True)
# allow the user to call setup
from . import setup  # noqa

# only import all other functionality if setup was successful
if _INSTANCE_SETUP:
    from lnschema_core import (  # noqa
        Featureset,
        File,
        Folder,
        Project,
        Run,
        Storage,
        Transform,
        User,
    )

    from . import _check_versions  # executes checks during import
    from . import dev  # noqa
    from . import schema  # noqa
    from . import types  # noqa
    from ._context import context  # noqa

    track = context._track  # noqa
    from lamin_logger import logger as _logger

    # this needs to follow on the import right now
    _logger.success(f"Loaded instance: {_lamindb_setup.settings.instance.identifier}")
    _logger.hint(f"Running lamindb {__version__}")

    from . import _amend_file  # noqa
    from . import _amend_folder  # noqa
    from ._delete import delete  # noqa
    from ._parse import parse  # noqa
    from ._save import save  # noqa
    from ._select import select  # noqa
    from ._settings import settings
    from ._view import view  # noqa

    add = save  # backward compat
