from typing import List

from lnschema_core import File, Folder

from .dev.db._select import select


# exposed to users as Folder.subset()
def subset(self: Folder, *, prefix: str, **fields) -> List[File]:
    """Get files via relative path to folder."""
    # ensure is actual folder, not virtual one
    if self.key is None:
        raise ValueError(
            ".get() is only defined for real folders, not virtual ones"
            "you can access files via .files or by refining with queries"
        )
    files = (
        select(File, **fields).where(File.key.startswith(self.key + "/" + prefix)).all()
    )
    return files


Folder.subset = subset
