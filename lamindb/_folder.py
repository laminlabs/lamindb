from itertools import islice
from pathlib import Path, PurePath
from typing import List, Optional, Union

import lndb
from lamin_logger import logger
from lndb_storage import UPath
from lnschema_core import File, Folder, Run
from lnschema_core.dev._storage import filepath_from_file_or_folder

from ._file import get_relative_path_to_directory
from .dev.db._select import select


def get_folder_kwargs_from_data(
    folder: Union[Path, UPath, str],
    *,
    name: Optional[str] = None,
    key: Optional[str] = None,
    source: Optional[Run] = None,
):
    folderpath = UPath(folder)
    cloudpath = folderpath if isinstance(folderpath, UPath) else None
    localpath = None if isinstance(folderpath, UPath) else folderpath
    folder_privates = dict(
        _local_filepath=localpath,
        _cloud_filepath=cloudpath,
    )
    if key is None:
        key = folderpath.name.rstrip("/")
        logger.hint(f"using key = {key}")

    # TODO: UPath doesn't list the first level files and dirs with "*"
    pattern = "" if isinstance(folderpath, UPath) else "*"

    files = []
    for filepath in folderpath.rglob(pattern):
        if filepath.is_file():
            relpath = get_relative_path_to_directory(filepath, folderpath)
            filekey = folderpath.name.rstrip("/") + "/" + relpath.as_posix()
            file = File(filepath, source=source, key=filekey)
            files.append(file)

    folder_kwargs = dict(
        name=folderpath.name if name is None else name,
        key=key,
        storage_id=lndb.settings.storage.id,
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
    folder: Folder, relpath: Union[str, Path, List[Union[str, Path]]], **fields
):
    """Get files via relative path to folder."""
    # ensure is actual folder, not virtual one
    if folder.key is None:
        raise ValueError(
            ".get() is only defined for real folders, not virtual ones"
            "you can access files via .files or by refining with queries"
        )

    if isinstance(relpath, List):
        relpaths = [PurePath(i) for i in relpath]
    else:
        abspath = relpath_to_abspath(folder=folder, relpath=PurePath(relpath))
        if abspath.is_dir():
            abspaths_files = list_files_in_directory(abspath)
            relpaths = [
                get_relative_path_to_directory(filepath=path, directory=folder.key)
                for path in abspaths_files
            ]
        else:
            relpaths = [PurePath(relpath)]

    file_keys = [relative_path_to_key(folder=folder, relpath=i) for i in relpaths]

    files = select(File, **fields).where(File.key.in_(file_keys)).all()
    return files


def list_files_in_directory(dirpath: Union[Path, UPath]):
    """List all files recursively from a directory."""
    return [path for path in dirpath.rglob("*") if path.is_file()]


def relative_path_to_key(folder: Folder, relpath: PurePath):
    """Convert a relative path of folder to an absolute path."""
    key = (PurePath(folder.key) / relpath.parent / relpath.name).as_posix()
    return key


def relpath_to_abspath(folder: Folder, relpath: PurePath):
    return filepath_from_file_or_folder(folder) / relpath
