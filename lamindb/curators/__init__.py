"""Curators.

.. autosummary::
   :toctree: .

   DataFrameCurator
   AnnDataCurator
   MuDataCurator
   SpatialDataCurator
   TiledbsomaExperimentCurator

Modules.

.. autosummary::
   :toctree: .

   core

"""

from ._legacy import (  # backward compat
    CellxGeneAnnDataCatManager,
    PertAnnDataCatManager,
)
from .core import (
    AnnDataCurator,
    DataFrameCurator,
    MuDataCurator,
    SpatialDataCurator,
    TiledbsomaExperimentCurator,
)

__all__ = [
    "CellxGeneAnnDataCatManager",
    "PertAnnDataCatManager",
    "AnnDataCurator",
    "DataFrameCurator",
    "MuDataCurator",
    "SpatialDataCurator",
    "TiledbsomaExperimentCurator",
]
