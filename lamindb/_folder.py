from itertools import islice
from pathlib import Path, PurePath
from typing import List, Optional, Union

from lndb.dev import UPath
from lnschema_core import DFolder as lns_DFolder
from lnschema_core import DObject as lns_DObject
from lnschema_core import Run, Storage

from ._record import get_storage_root_and_root_str
from .dev._core import filepath_from_dfolder, get_name_suffix_from_filepath  # noqa
from .dev.db._add import write_objectkey
from .dev.db._select import select


def get_dfolder_kwargs_from_data(
    folder: Union[Path, UPath, str],
    *,
    name: Optional[str] = None,
    source: Optional[Run] = None,
):
    folderpath = UPath(folder)
    cloudpath = folderpath if isinstance(folderpath, UPath) else None
    localpath = None if isinstance(folderpath, UPath) else folderpath
    dfolder_privates = dict(
        _local_filepath=localpath,
        _cloud_filepath=cloudpath,
    )

    dobjects = []
    for f in folderpath.glob("**/*"):
        if f.is_file():
            dobj = lns_DObject(f, source=source)
            write_objectkey(dobj)
            dobjects.append(dobj)

    dfolder_kwargs = dict(
        name=folderpath.name if name is None else name,
        dobjects=dobjects,
    )
    return dfolder_kwargs, dfolder_privates


# Exposed to users as DFolder.tree()
def tree(
    dir_path: Union[Path, UPath, str],
    level: int = -1,
    limit_to_directories: bool = False,
    length_limit: int = 1000,
):
    """Given a directory Path object print a visual tree structure.

    Adapted from: https://stackoverflow.com/questions/9727673/list-directory-tree-structure-in-python  # noqa
    """
    space = "    "
    branch = "│   "
    tee = "├── "
    last = "└── "

    dir_path = UPath(dir_path)
    files = 0
    directories = 0

    def inner(dir_path: Path, prefix: str = "", level=-1):
        nonlocal files, directories
        if not level:
            return  # 0, stop iterating
        if limit_to_directories:
            contents = [d for d in dir_path.iterdir() if d.is_dir()]
        else:
            contents = list(dir_path.iterdir())
        pointers = [tee] * (len(contents) - 1) + [last]
        for pointer, path in zip(pointers, contents):
            if path.is_dir():
                yield prefix + pointer + path.name
                directories += 1
                extension = branch if pointer == tee else space
                yield from inner(path, prefix=prefix + extension, level=level - 1)
            elif not limit_to_directories:
                yield prefix + pointer + path.name
                files += 1

    print(dir_path.name)
    iterator = inner(dir_path, level=level)
    for line in islice(iterator, length_limit):
        print(line)
    if next(iterator, None):
        print(f"... length_limit, {length_limit}, reached, counted:")
    print(f"\n{directories} directories" + (f", {files} files" if files else ""))


# Exposed to users as DFolder.get()
def get_dobject(
    dfolder: lns_DFolder, relpath: Union[str, Path, List[Union[str, Path]]], **fields
):
    """Get dobjects via relative path to dfolder."""
    if isinstance(relpath, List):
        relpaths = [PurePath(i) for i in relpath]
    else:
        relpaths = [PurePath(relpath)]

    dobject_objectkeys = get_dobject_objectkey(dfolder=dfolder, relpaths=relpaths)
    return select_by_objectkey(dobject_objectkeys=dobject_objectkeys, **fields)


def get_dobject_objectkey(dfolder: lns_DFolder, relpaths: List[PurePath]):
    """Get dobject via relative path to dfolder."""
    keys = []
    for rpath in relpaths:
        name, _ = get_name_suffix_from_filepath(rpath)
        keys.append(str(PurePath(dfolder._objectkey) / rpath.parent / name))

    # check relative path exists (slow for cloud path)
    # if not dobject_path.exists():
    #     raise AssertionError(f"Unable to locate dobject at {dobject_path}")
    # dobject_path = filepath_from_dfolder(dfolder) / relpath

    return keys


def select_by_objectkey(dobject_objectkeys: List[str], **fields):
    # query for unique comb of (_dobjectkey, storage, suffix)
    dobjects = (
        select(lns_DObject, **fields)
        .where(lns_DObject._objectkey.in_(dobject_objectkeys))
        .join(Storage, root=get_storage_root_and_root_str()[1])
        .all()
    )
    # TODO: return dobjects in the same order as the dobject_objectkeys
    return dobjects
