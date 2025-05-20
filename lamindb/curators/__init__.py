"""Curators.

.. autosummary::
   :toctree: .

   DataFrameCurator
   AnnDataCurator
   SpatialDataCurator
   MuDataCurator

Modules.

.. autosummary::
   :toctree: .

   core

"""

from ._legacy import (  # backward compat
    CellxGeneAnnDataCatManager,
    PertAnnDataCatManager,
)
from .core import AnnDataCurator, DataFrameCurator, MuDataCurator, SpatialDataCurator
