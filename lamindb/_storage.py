from lamindb_setup.dev.upath import UPath
from lnschema_core import Storage

from lamindb.dev.storage.file import _str_to_path


def root_as_path(self) -> UPath:
    return _str_to_path(self.root)


setattr(Storage, "root_as_path", root_as_path)
