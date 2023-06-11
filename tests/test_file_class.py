import shutil
from pathlib import Path

import pytest
from lnschema_core.models import File

# how do we properly abstract out the default storage variable?
# currently, we're only mocking it through `default_storage` as
# set in conftest.py


@pytest.fixture(scope="module")
def create_dummy_files():
    # create a file within the default storage
    dir1 = Path("./default_storage") / Path("my_project")
    dir1.mkdir(parents=True)
    filepath1 = dir1 / "my_file1.csv"
    filepath1.write_text("a")
    # create a file outside the default storage
    dir2 = Path("./outside_storage") / Path("my_project")
    dir2.mkdir(parents=True)
    filepath2 = dir2 / "my_file2.csv"
    filepath2.write_text("b")
    yield
    shutil.rmtree(dir1)
    shutil.rmtree(dir2)


def test_init_from_filepath(create_dummy_files):
    # file within storage, key=None, name=None
    file = File("./default_storage/my_project/my_file1.csv")
    assert file.name == "my_file1.csv"
    assert file.suffix == ".csv"
    assert file.key == "my_project/my_file1.csv"
    assert file.storage.root == Path("./default_storage").resolve().as_posix()
    assert file.hash == "DMF1ucDxtqgxw5niaXcmYQ"
    # file outside storage, key=None, name=None
    file = File("./outside_storage/my_project/my_file2.csv")
    assert file.name == "my_file2.csv"
    assert file.suffix == ".csv"
    assert file.key is None
    assert file.storage.root == Path("./default_storage").resolve().as_posix()
    assert file.hash == "kutf_uauL-w61xx3dTFXjw"
    # file within storage, key=None, name!=None
    file = File("./default_storage/my_project/my_file1.csv", name="hello")
    assert file.name == "hello"
    assert file.suffix == ".csv"
    assert file.key == "my_project/my_file1.csv"
    assert file.storage.root == Path("./default_storage").resolve().as_posix()
    assert file.hash == "DMF1ucDxtqgxw5niaXcmYQ"
    # file outside storage, key=None, name!=None
    file = File("./outside_storage/my_project/my_file2.csv", name="hello")
    assert file.name == "hello"
    assert file.suffix == ".csv"
    assert file.key is None
    assert file.storage.root == Path("./default_storage").resolve().as_posix()
    assert file.hash == "kutf_uauL-w61xx3dTFXjw"
