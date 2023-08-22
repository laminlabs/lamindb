from lamindb_setup import settings
from lamindb_setup.dev.upath import UPath
from lnschema_core import Storage


def root_as_path(self) -> UPath:
    return settings.storage.to_path(self.root)


setattr(Storage, "root_as_path", root_as_path)
