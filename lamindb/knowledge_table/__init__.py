"""Knowledge tables.

Feature tables from bionty:

.. autosummary::
   :toctree: .

   Gene
   Protein
   CellMarker

Lookup knowledge table identifiers:

.. autosummary::
   :toctree: .

   lookup
"""


from ._core import CellMarker, Gene, Protein
from ._lookup import lookup
