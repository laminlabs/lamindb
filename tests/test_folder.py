import pytest

import lamindb as ln
from lamindb.dev.storage import UPath


def test_folder_tree():
    # virtual folder
    new_dir = UPath("./test-folder")
    new_dir.mkdir(exist_ok=True)
    with open(new_dir / "file-in-folder.txt", "w") as f:
        f.write("file-in-folder")
    with pytest.raises(NotImplementedError):
        folder = ln.Folder(new_dir)
        folder.tree()

    # real folder
    new_dir = UPath(f"{ln.setup.settings.storage.root}/folder")
    new_dir.mkdir(exist_ok=True)
    with open(new_dir / "file-in-folder.txt", "w") as f:
        f.write("file-in-folder")
    folder = ln.Folder(new_dir)
    folder.tree()
