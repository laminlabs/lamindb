import os
from collections import defaultdict
from typing import Iterable

from lamindb_setup import settings as setup_settings
from lnschema_core.models import File, Storage


def view_tree(
    cls,
    level: int = -1,
    limit_to_directories: bool = False,
    length_limit: int = 1000,
    max_files_per_dir_per_type: int = 7,
) -> None:
    """{}"""
    if cls.__class__.__name__ == "QuerySet":
        print("queryset")
        qs = cls
        storage_ids = qs.list("storage_id")
    elif cls == File:
        print("file")
        qs = cls.filter(storage_id=setup_settings.storage.id).all()
        storage_ids = Storage.filter().list("id")
    else:
        print("else")
        return
    storages = Storage.filter().all()
    storage_roots = {
        storage_id: storages.get(id=storage_id).root for storage_id in storage_ids
    }
    keys = set()
    for file in qs:
        root = storage_roots.get(file.storage_id, "")
        keys.add(f"{root}/{file.key}")

    _view_tree(
        keys=keys,
        level=level,
        only_dirs=limit_to_directories,
        limit=length_limit,
        max_files_per_dir_per_type=max_files_per_dir_per_type,
    )


def _view_tree(
    keys: Iterable[str],
    *,
    level: int = -1,
    only_dirs: bool = False,
    limit: int = 1000,
    max_files_per_dir_per_type: int = 7,
) -> None:
    # Create a nested dictionary from keys
    def tree():
        return defaultdict(tree)

    root = tree()

    n_files = 0
    n_directories = 0
    suffixes = set()

    for key in keys:
        parts = key.split("/")
        node = root
        for part in parts:
            node = node[part]
            if node == {}:
                n_files += 1
                suffix = os.path.splitext(part)[1]
                if suffix:
                    suffixes.add(suffix)
            else:
                n_directories += 1

    # Function to print the tree
    def print_tree(node, prefix="", depth=0, count=[0], n_files_per_dir_per_type=None):
        if n_files_per_dir_per_type is None:
            n_files_per_dir_per_type = defaultdict(int)

        if level != -1 and depth > level:
            return
        for name, child in node.items():
            if count[0] >= limit:
                return
            if only_dirs and child == {}:
                continue
            suffix = os.path.splitext(name)[1]
            n_files_per_dir_per_type[suffix] += 1
            if (
                depth > 0
                and n_files_per_dir_per_type[suffix] > max_files_per_dir_per_type
            ):
                continue
            new_prefix = prefix + ("├── " if name != list(node.keys())[-1] else "└── ")
            print(new_prefix + name)
            count[0] += 1
            if child:
                print_tree(
                    child,
                    prefix + ("│   " if name != list(node.keys())[-1] else "    "),
                    depth + 1,
                    count,
                    (
                        defaultdict(int) if depth == 0 else n_files_per_dir_per_type
                    ),  # Reset the counter for each directory
                )

    suffix_message = f" with suffixes {', '.join(suffixes)}" if n_files > 0 else ""
    print(f"{n_directories} directories, {n_files} files{suffix_message}")
    print_tree(root)
