"""Integrations.

Modules
-------

.. autosummary::
   :toctree: .

   lightning

Classes
-------

.. autoclass:: Reader

Functions
---------

.. autofunction:: save_vitessce_config
.. autofunction:: save_tiledbsoma_experiment
.. autofunction:: curate_from_croissant
.. autofunction:: import_db
.. autofunction:: link

"""

from ._croissant import curate_from_croissant
from ._vitessce import save_vitessce_config

__all__ = [
    "lightning",
    "Reader",
    "save_tiledbsoma_experiment",
    "curate_from_croissant",
    "save_vitessce_config",
    "import_db",
    "link",
    "upsert",
    "materialize",
]

_NOTION_EXPORTS = {"Reader", "import_db", "link", "upsert", "materialize"}


def __getattr__(name: str):
    """Lazy-import heavy symbols to avoid loading storage/lamindb at package import."""
    if name == "save_tiledbsoma_experiment":
        from lamindb.core.storage import save_tiledbsoma_experiment

        return save_tiledbsoma_experiment

    if name in _NOTION_EXPORTS:
        from . import notion

        return getattr(notion, name)

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
