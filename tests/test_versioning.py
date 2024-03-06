import lamindb as ln
import pandas as pd
import pytest
from lamindb import UPath
from lamindb.core.versioning import get_new_path_from_uid, set_version


@pytest.fixture(scope="module")
def df1():
    return pd.DataFrame({"feat1": [1, 2]})


@pytest.fixture(scope="module")
def df2():
    return pd.DataFrame({"feat1": [2, 3]})


def test_set_version():
    # all remaining lines are covered in notebooks
    with pytest.raises(ValueError):
        set_version(None, "1.2")
    assert set_version(None, "0") == "1"
    assert set_version(None, "1") == "2"
    assert set_version("1.2.3", "0") == "1.2.3"
    assert set_version("1.2.3") == "1.2.3"


def test_add_to_version_family(df1, df2):
    artifact1 = ln.Artifact.from_df(df1, description="test1")
    artifact1.save()
    artifact2 = ln.Artifact.from_df(df2, description="test2")
    artifact2.save()
    assert (
        artifact1.uid[: artifact1._len_stem_uid]
        != artifact2.uid[: artifact2._len_stem_uid]
    )
    artifact2.add_to_version_family(artifact1)
    assert (
        artifact1.uid[: artifact1._len_stem_uid]
        == artifact2.uid[: artifact2._len_stem_uid]
    )
    assert (
        artifact1.path.name[: artifact1._len_stem_uid]
        == artifact2.path.name[: artifact2._len_stem_uid]
    )
    artifact1.delete(permanent=True)
    artifact2.delete(permanent=True)


def test_get_new_path_from_uid():
    # test cloud path as it has different behavior than local path
    with open("test_new_path.txt", "w") as f:
        f.write("test_new_path")
    old_path = UPath("s3://lamindata/.lamindb/test_new_path.txt")
    old_path.upload_from("./test_new_path.txt")
    assert old_path.exists()
    new_path = get_new_path_from_uid(
        old_path=old_path,
        old_uid="test_new_path",
        new_uid="test_new_path2",
    )
    assert new_path == "test_new_path2.txt"
    new_path = old_path.rename(new_path)
    assert new_path.exists()
    assert str(new_path) == "s3://lamindata/.lamindb/test_new_path2.txt"
    assert not old_path.exists()
    new_path.unlink()
    UPath("./test_new_path.txt").unlink()
