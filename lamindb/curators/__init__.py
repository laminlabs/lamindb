"""Curators.

.. autosummary::
   :toctree: .

   DataFrameCurator
   AnnDataCurator
   MuDataCurator
   SpatialDataCurator

Modules.

.. autosummary::
   :toctree: .

   core

"""

from ._legacy import PertAnnDataCatManager  # backward compat
from .core import AnnDataCurator, DataFrameCurator, MuDataCurator, SpatialDataCurator
