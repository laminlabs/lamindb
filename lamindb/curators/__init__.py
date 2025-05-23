"""Curators.

.. autosummary::
   :toctree: .

   DataFrameCurator
   AnnDataCurator
   MuDataCurator
   SpatialDataCurator
   SomaExperimentCurator

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
    SomaExperimentCurator,
    SpatialDataCurator,
)

__all__ = [
    "CellxGeneAnnDataCatManager",
    "PertAnnDataCatManager",
    "AnnDataCurator",
    "DataFrameCurator",
    "MuDataCurator",
    "SpatialDataCurator",
    "SomaExperimentCurator",
]
