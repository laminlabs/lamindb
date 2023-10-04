from lamindb_setup.dev._docs import doc_args
from lamindb_setup.dev.upath import UPath, create_path
from lnschema_core import Storage


def root_as_path(self) -> UPath:
    return create_path(self.root)


@property  # type: ignore
@doc_args(Storage.path.__doc__)
def path(self) -> UPath:
    """{}"""
    return create_path(self.root)


setattr(Storage, "root_as_path", root_as_path)
setattr(Storage, "path", path)
