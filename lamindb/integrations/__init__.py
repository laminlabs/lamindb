"""Integrations.

.. autosummary::
   :toctree: .

   save_vitessce_config
   save_tiledbsoma_experiment
   curate_from_croissant
   lightning
"""

from typing import Any


def __getattr__(attr_name: str) -> Any:
    # Defers import until accessed to avoid requiring PyTorch Lightning
    if attr_name == "lightning":
        from lamindb.integrations import lightning

        return lightning
    raise AttributeError(f"module has no attribute {attr_name!r}")


from lamindb.core.storage import save_tiledbsoma_experiment

from ._croissant import curate_from_croissant
from ._vitessce import save_vitessce_config


def __dir__():
    # Makes lazy imports discoverable to dir() to enable autocomplete including lazy modules
    return __all__


__all__ = [
    "lightning",
    "save_tiledbsoma_experiment",
    "curate_from_croissant",
    "save_vitessce_config",
]
