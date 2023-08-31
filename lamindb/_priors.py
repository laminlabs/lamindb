from functools import cached_property

from lnschema_core import Modality, User

Modality(name="meta").save()


# might be a bad design idea & dispear again
class Priors:
    @cached_property
    def modalities(self):
        return Modality.lookup()

    @cached_property
    def users(self):
        return User.lookup()


priors = Priors()
