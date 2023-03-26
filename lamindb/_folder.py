from itertools import islice
from pathlib import Path, PurePath
from typing import List, Optional, Union

from lndb_storage import UPath
from lnschema_core import File as lns_File
from lnschema_core import Folder as lns_Folder
from lnschema_core import Run, Storage

from ._record import get_storage_root_and_root_str
from .dev._core import filepath_from_folder, get_name_suffix_from_filepath
from .dev.db._add import filepath_to_relpath, write_objectkey
from .dev.db._select import select


def get_folder_kwargs_from_data(
    folder: Union[Path, UPath, str],
    *,
    name: Optional[str] = None,
    source: Optional[Run] = None,
):
    folderpath = UPath(folder)
    cloudpath = folderpath if isinstance(folderpath, UPath) else None
    localpath = None if isinstance(folderpath, UPath) else folderpath
    folder_privates = dict(
        _local_filepath=localpath,
        _cloud_filepath=cloudpath,
    )

    # TODO: UPath doesn't list the first level files and dirs with "*"
    pattern = "" if isinstance(folderpath, UPath) else "*"

    files = []
    for f in folderpath.rglob(pattern):
        if f.is_file():
            dobj = lns_File(f, source=source)
            write_objectkey(dobj)
            files.append(dobj)

    folder_kwargs = dict(
        name=folderpath.name if name is None else name,
        files=files,
    )
    return folder_kwargs, folder_privates


# Exposed to users as Folder.tree()
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
        # this is needed so that the passed folder is not listed
        contents = [
            i
            for i in dir_path.iterdir()
            if i.as_posix().rstrip("/") != dir_path.as_posix().rstrip("/")
        ]
        if limit_to_directories:
            contents = [d for d in contents if d.is_dir()]
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


# Exposed to users as Folder.get()
def get_file(
    folder: lns_Folder, relpath: Union[str, Path, List[Union[str, Path]]], **fields
):
    """Get files via relative path to folder."""
    if isinstance(relpath, List):
        relpaths = [PurePath(i) for i in relpath]
    else:
        abspath = relpath_to_abspath(folder=folder, relpath=PurePath(relpath))
        if abspath.is_dir():
            abspaths_files = list_files_from_dir(abspath)
            root, root_str = get_storage_root_and_root_str(filepath_from_folder(folder))
            relpaths = [
                filepath_to_relpath(root=root, root_str=root_str, filepath=i)
                for i in abspaths_files
            ]
        else:
            relpaths = [PurePath(relpath)]

    file_objectkeys = [relpath_to_objectkey(folder=folder, relpath=i) for i in relpaths]

    return select_by_objectkey(file_objectkeys=file_objectkeys, **fields)


def list_files_from_dir(dirpath: Union[Path, UPath]):
    """List all files recursively from a directory."""
    return [i for i in dirpath.rglob("*") if i.is_file()]


def relpath_to_objectkey(folder: lns_Folder, relpath: PurePath):
    """Convert a relative path of folder to an absolute path."""
    name, _ = get_name_suffix_from_filepath(relpath)
    objectkey = str(PurePath(folder._objectkey) / relpath.parent / name)
    return objectkey


def relpath_to_abspath(folder: lns_Folder, relpath: PurePath):
    return filepath_from_folder(folder) / relpath


def select_by_objectkey(file_objectkeys: List[str], **fields):
    # query for unique comb of (_filekey, storage, suffix)
    files = (
        select(lns_File, **fields)
        .where(lns_File._objectkey.in_(file_objectkeys))
        .join(Storage, root=get_storage_root_and_root_str()[1])
        .all()
    )
    # TODO: return files in the same order as the file_objectkeys
    return files
