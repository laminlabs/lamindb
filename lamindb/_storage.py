from lamindb_setup.dev.upath import UPath
from lnschema_core import Storage


def root_as_path(self) -> UPath:
    return UPath(self.root)


setattr(Storage, "root_as_path", root_as_path)
