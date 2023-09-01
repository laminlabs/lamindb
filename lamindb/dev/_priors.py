from functools import cached_property

from lnschema_core import Modality

from .._save import save
from ._settings import settings

MODALITIES = ["meta"]
NUMBER_TYPE = "number"


# might be a bad design idea & disappear again
class Priors:
    @cached_property
    def modalities(self):
        all_defaults = Modality.filter(name__in=MODALITIES).all()
        # update the registry
        if len(all_defaults) != len(MODALITIES):
            upon_create_search_names = settings.upon_create_search_names
            settings.upon_create_search_names = False
            modalities = [Modality(name=name) for name in MODALITIES]
            settings.upon_create_search_names = upon_create_search_names
            save(modalities)
        return Modality.lookup()


priors = Priors()
