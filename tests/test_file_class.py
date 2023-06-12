import shutil
from pathlib import Path

import pytest
from lnschema_core.models import File

# how do we properly abstract out the default storage variable?
# currently, we're only mocking it through `default_storage` as
# set in conftest.py


@pytest.fixture(
    scope="module", params=[(True, "./default_storage/"), (False, "./outside_storage/")]
)
def get_test_filepaths(request):
    isin_default_storage: bool = request.param[0]
    root_dir: Path = Path(request.param[1])
    test_folder = root_dir / "my_folder/"
    test_folder.mkdir(parents=True)
    test_filepath = test_folder / "my_file.csv"
    test_filepath.write_text("a")
    # return a boolean indicating whether test filepath is in default storage
    # and the test filepath
    yield (isin_default_storage, test_filepath)
    shutil.rmtree(test_folder)


# this tests the basic (non-provenance-related) metadata
@pytest.mark.parametrize("name", [None, "my name"])
def test_init_from_filepath_basic_fields(get_test_filepaths, name):
    isin_default_storage = get_test_filepaths[0]
    test_filepath = get_test_filepaths[1]
    file = File(test_filepath, name=name)
    assert file.name == test_filepath.name if name is None else file.name == name
    assert file.suffix == ".csv"
    assert (
        file.key == "my_folder/my_file.csv"
        if isin_default_storage
        else file.key is None
    )
    assert file.storage.root == Path("./default_storage").resolve().as_posix()
    assert file.hash == "DMF1ucDxtqgxw5niaXcmYQ"
