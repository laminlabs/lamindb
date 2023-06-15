import shutil
from pathlib import Path

import pytest

import lamindb as ln


@pytest.fixture(
    scope="module", params=[(True, "./default_storage/"), (False, "./outside_storage/")]
)
def get_test_folderpaths(request):
    isin_default_storage: bool = request.param[0]
    root_dir: Path = Path(request.param[1])
    test_folder = root_dir / "my_folder/"
    test_folder.mkdir(parents=True)
    test_filepath = test_folder / "my_file.csv"
    test_filepath.write_text("a")
    # return a boolean indicating whether test filepath is in default storage
    # and the test filepath
    yield (isin_default_storage, test_folder)
    shutil.rmtree(test_folder)


def test_folder_init_from_path(get_test_folderpaths):
    isin_default_storage = get_test_folderpaths[0]
    test_folderpath = get_test_folderpaths[1]
    if not isin_default_storage:
        with pytest.raises(RuntimeError):
            ln.Folder(test_folderpath)
    else:
        folder = ln.Folder(test_folderpath)
        # also execute tree
        folder.tree()
