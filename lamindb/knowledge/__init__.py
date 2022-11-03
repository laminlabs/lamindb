"""Knowledge.

Feature tables from bionty:

.. autosummary::
   :toctree: .

   Species
   Gene
   Protein
   CellMarker

Lookup knowledge table identifiers:

.. autosummary::
   :toctree: .

   lookup
"""


from ._core import CellMarker, Gene, Protein, Species
from ._lookup import lookup
