from itertools import islice
from pathlib import Path
from typing import Optional, Union

import lamindb_setup
from lamin_logger import logger
from lnschema_core import ids
from lnschema_core.models import File, Folder, Run

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
        raise RuntimeError(
            "Only folders in the default storage can be registered!\n"
            "You can either move your folder into the current default storage"
            "or add a new default storage through ln.setup.set.storage()"
        )
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
    folder: Folder,
    level: int = -1,
    limit_to_directories: bool = False,
    length_limit: int = 1000,
):
    """Given a directory Path object print a visual tree structure.

    Adapted from: https://stackoverflow.com/questions/9727673/list-directory-tree-structure-in-python  # noqa
    """
    if folder.key is None:
        raise RuntimeError("Virtual folders do not have a tree structure")

    space = "    "
    branch = "│   "
    tee = "├── "
    last = "└── "

    dir_path = UPath(folder.path())
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

    folder_tree = f"{dir_path.name}"
    iterator = inner(dir_path, level=level)
    for line in islice(iterator, length_limit):
        folder_tree += f"\n{line}"
    if next(iterator, None):
        folder_tree += f"... length_limit, {length_limit}, reached, counted:"
    print(folder_tree)
    print(f"\n{directories} directories" + (f", {files} files" if files else ""))


def init_folder(folder: Folder, *args, **kwargs):
    # Below checks for the Django-internal call in from_db()
    # it'd be better if we could avoid this, but not being able to create a File
    # from data with the default constructor renders the central class of the API
    # essentially useless
    # The danger below is not that a user might pass as many args (12 of it), but rather
    # that at some point the Django API might change; on the other hand, this
    # condition of for calling the constructor based on kwargs should always
    # stay robust
    if len(args) == len(folder._meta.concrete_fields):
        super(Folder, folder).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) != 1 and "files" not in kwargs:
        raise ValueError("Either provide path as arg or provide files as kwarg!")
    if len(args) == 1:
        path: Optional[Union[Path, UPath, str]] = args[0]
    else:
        path = None
    name: Optional[str] = kwargs.pop("name") if "name" in kwargs else None
    key: Optional[str] = kwargs.pop("key") if "key" in kwargs else None
    files: Optional[str] = kwargs.pop("files") if "files" in kwargs else None
    if len(kwargs) != 0:
        raise ValueError(f"These kwargs are not permitted: {kwargs}")

    if path is not None:
        kwargs, privates = get_folder_kwargs_from_data(
            path=path,
            name=name,
            key=key,
        )
        files = kwargs.pop("files")
    else:
        kwargs = dict(name=name)
    kwargs["id"] = ids.base62_20()
    super(Folder, folder).__init__(**kwargs)
    if path is not None:
        folder._local_filepath = privates["local_filepath"]
        folder._cloud_filepath = privates["cloud_filepath"]
        folder._files = files
