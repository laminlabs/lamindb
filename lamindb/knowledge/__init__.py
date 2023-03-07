"""Knowledge.

Feature tables from bionty:

.. autosummary::
   :toctree: .

   Species
   Gene
   Protein
   CellMarker
   Tissue
   CellType
   Disease

Lookup knowledge table identifiers:

.. autosummary::
   :toctree: .

   lookup
"""


from lamin_logger import logger

from ._core import CellMarker, CellType, Disease, Gene, Protein, Species, Tissue
from ._lookup import lookup

logger.warning(
    "The lamindb.knowledge API is deprecated.\nPlease, use bionty directly.\nYou can"
    " replace all occurances of ln.knowledge with bionty without breaking changes!"
)
