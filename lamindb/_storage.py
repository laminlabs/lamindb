from lamindb_setup.dev.upath import UPath, create_path
from lnschema_core import Storage


def root_as_path(self) -> UPath:
    return create_path(self.root)


setattr(Storage, "root_as_path", root_as_path)
