"""Integrations.

.. autosummary::
   :toctree: .

   save_vitessce_config
   save_tiledbsoma_experiment
   curate_from_croissant
   LightningCallback
"""

from typing import Any


def __getattr__(attr_name: str) -> Any:
    if attr_name == "LightningCallback":
        from ._lightning import LightningCallback

        return LightningCallback


from lamindb.core.storage import save_tiledbsoma_experiment

from ._croissant import curate_from_croissant
from ._vitessce import save_vitessce_config

__all__ = [
    "LightningCallback",
    "save_tiledbsoma_experiment",
    "curate_from_croissant",
    "save_vitessce_config",
]
