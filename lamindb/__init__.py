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

Tracking data by features:

.. autosummary::
   :toctree: .

   Features

Query & manipulate data:

.. autosummary::
   :toctree: .

   select
   add
   delete

Manipulate data with open session:

.. autosummary::
   :toctree: .

   Session

Utility functions:

.. autosummary::
   :toctree: .

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
   dev
"""

__version__ = "0.39rc1"  # denote a release candidate for 0.1.0 with 0.1rc1

# prints warning of python versions
from lamin_logger import py_version_warning as _py_version_warning
from lndb._check_instance_setup import check_instance_setup as _check_instance_setup

_py_version_warning("3.8", "3.10")

_INSTANCE_SETUP = _check_instance_setup(from_lamindb=True)

# allow the user to call setup
from . import setup  # noqa

# only import all other functionality if setup was successful
if _INSTANCE_SETUP:
    from lndb_storage import subset
    from lndb_storage.object import lazy
    from lnschema_core import (  # noqa
        Features,
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
    from ._context import context  # noqa

    track = context._track  # noqa
    from . import _amend_file  # noqa
    from . import _amend_folder  # noqa
    from ._delete import delete  # noqa
    from ._nb import nb  # noqa
    from ._settings import settings
    from ._view import view  # noqa
    from .dev.db import Session  # noqa
    from .dev.db._add import add  # noqa
    from .dev.db._select import select  # noqa
