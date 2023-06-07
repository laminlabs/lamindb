from itertools import islice
from pathlib import Path
from typing import Optional, Union

import lamindb_setup
from lamin_logger import logger
from lnschema_core import File, Run
from lnschema_core import ids as id_generator

from lamindb.dev.storage import UPath

from ._file import (
    get_check_path_in_storage,
    get_relative_path_to_directory,
    get_relative_path_to_root,
)


def get_folder_kwargs_from_data(
    path: Union[Path, UPath, str],
    *,
    name: Optional[str] = None,
    key: Optional[str] = None,
    run: Optional[Run] = None,
):
    folderpath = UPath(path)
    check_path_in_storage = get_check_path_in_storage(folderpath)

    if key is None and check_path_in_storage:
        folder_key = get_relative_path_to_root(path=folderpath).as_posix()
    elif key is None:
        folder_key = id_generator.folder()
    else:
        folder_key = key

    # always sanitize by stripping a trailing slash
    folder_key = folder_key.rstrip("/")

    logger.hint(f"using storage key = {folder_key}")

    # TODO: UPath doesn't list the first level files and dirs with "*"
    pattern = "" if isinstance(folderpath, UPath) else "*"

    files = []
    for filepath in folderpath.rglob(pattern):
        if filepath.is_file():
            relative_path = get_relative_path_to_directory(filepath, folderpath)
            file_key = folder_key + "/" + relative_path.as_posix()
            files.append(File(filepath, run=run, key=file_key))

    logger.hint(f"-> n_files = {len(files)}")

    kwargs = dict(
        name=folderpath.name if name is None else name,
        key=folder_key,
        storage_id=lamindb_setup.settings.storage.id,
        files=files,
    )
    privates = dict(
        local_filepath=folderpath if isinstance(folderpath, UPath) else None,
        cloud_filepath=None if isinstance(folderpath, UPath) else folderpath,
    )
    return kwargs, privates


# exposed to users as Folder.tree()
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

    try:
        folder_tree = f"{dir_path.name}"
        iterator = inner(dir_path, level=level)
        for line in islice(iterator, length_limit):
            folder_tree += f"\n{line}"
        if next(iterator, None):
            folder_tree += f"... length_limit, {length_limit}, reached, counted:"
        print(folder_tree)
        print(f"\n{directories} directories" + (f", {files} files" if files else ""))
    except FileNotFoundError:
        raise NotImplementedError(
            "Tree display only works for real folders: folders that exist in storage."
        )
