from typing import List

from lamindb_setup import _USE_DJANGO
from lnschema_core import File, Folder

from .dev.db._select import select

Folder.__doc__ = """Folders: collections of files.

Real vs. virtual folders:

- A real LaminDB `Folder` has a 1:1 correspondence to a folder on a file system
  or in object storage, and a `.key` that is not `None`.
- A virtual LaminDB `Folder` is a mere way of grouping files. A file can be
  linked to multiple virtual folders, but only to one real folder.
"""


# exposed to users as Folder.subset()
def subset(self: Folder, *, prefix: str, **fields) -> List[File]:
    """Get files via relative path to folder."""
    # ensure is actual folder, not virtual one
    if self.key is None:
        raise ValueError(
            ".get() is only defined for real folders, not virtual ones"
            "you can access files via .files or by refining with queries"
        )
    if _USE_DJANGO:
        files = (
            select(File, **fields)
            .filter(key__startswith=self.key + "/" + prefix)
            .list()
        )
    else:
        files = (
            select(File, **fields)
            .where(File.key.startswith(self.key + "/" + prefix))
            .all()
        )
    return files


Folder.subset = subset
