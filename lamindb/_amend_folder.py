from pathlib import Path
from typing import List, Union

from lnschema_core import File, Folder
from upath import UPath

from ._select import select

Folder.__doc__ = """Folders: collections of files.

Real vs. virtual folders:

- A real LaminDB `Folder` has a 1:1 correspondence to a folder on a file system
  or in object storage, and a `.key` that is not `None`.
- A virtual LaminDB `Folder` is a mere way of grouping files. A file can be
  linked to multiple virtual folders, but only to one real folder.
"""


# exposed to users as Folder.subset()
def subset(folder: Folder, *, prefix: str, **fields) -> List[File]:
    """Get files via relative path to folder."""
    # ensure is actual folder, not virtual one
    if folder.key is None:
        raise ValueError(
            ".get() is only defined for real folders, not virtual ones"
            "you can access files via .files or by refining with queries"
        )
    files = (
        select(File, **fields).filter(key__startswith=folder.key + "/" + prefix).list()
    )
    return files


def path(folder: Folder) -> Union[Path, UPath]:
    """Path on storage."""
    from lamindb._file_access import filepath_from_file_or_folder

    return filepath_from_file_or_folder(folder)


def tree(
    folder: Folder,
    level: int = -1,
    limit_to_directories: bool = False,
    length_limit: int = 1000,
) -> None:
    """Print a visual tree structure."""
    from lamindb._folder import tree

    return tree(
        folder,
        level=level,
        limit_to_directories=limit_to_directories,
        length_limit=length_limit,
    )


def __init__(folder: Folder, *args, **kwargs):
    from lamindb._folder import init_folder

    init_folder(folder, *args, **kwargs)


def save(folder: Folder, *args, **kwargs) -> None:
    """Save the folder."""
    # only has attr _files if freshly initialized
    if hasattr(folder, "_files"):
        for file in folder._files:
            file.save()
    super(Folder, folder).save(*args, **kwargs)
    if hasattr(folder, "_files"):
        folder.files.set(folder._files)


Folder.subset = subset
Folder.save = save
Folder.tree = tree
Folder.path = path
