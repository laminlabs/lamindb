import lamindb as ln
import pandas as pd
import pytest
from lamindb import UPath
from lamindb.models._is_versioned import (
    bump_version,
    get_new_path_from_uid,
    set_version,
)


@pytest.fixture(scope="module")
def df1():
    return pd.DataFrame({"feat1": [1, 2]})


@pytest.fixture(scope="module")
def df2():
    return pd.DataFrame({"feat1": [2, 3]})


def test_set_version():
    # all remaining lines are covered in notebooks
    with pytest.raises(ValueError):
        set_version(None, "weird-version")
    assert set_version(None, "1.2") == "2"
    assert set_version(None, "0") == "1"
    assert set_version(None, "1") == "2"
    assert set_version("1.2.3", "0") == "1.2.3"
    assert set_version("1.2.3") == "1.2.3"


def test_bump_version():
    current_version_major_only = "2"
    current_version_major_minor = "2.1"
    weird_version = "weird-version"
    with pytest.raises(ValueError):
        bump_version(weird_version)
    assert bump_version(weird_version, behavior="ignore") == "?"
    assert bump_version(current_version_major_only, bump_type="major") == "3"
    assert bump_version(current_version_major_only, bump_type="minor") == "2.1"
    assert bump_version(current_version_major_minor, bump_type="major") == "3"
    assert bump_version(current_version_major_minor, bump_type="minor") == "2.2"


def test__add_to_version_family(df1, df2):
    artifact1 = ln.Artifact.from_df(df1, description="test1")
    artifact1.save()
    artifact2 = ln.Artifact.from_df(df2, description="test2")
    artifact2.save()
    assert (
        artifact1.uid[: artifact1._len_stem_uid]
        != artifact2.uid[: artifact2._len_stem_uid]
    )
    artifact2._add_to_version_family(artifact1)
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


def test_latest_version_and_get():
    # build one version family
    transform_v1 = ln.Transform(description="Introduction").save()
    assert transform_v1.is_latest
    assert transform_v1.version is None
    # pass the latest version, also vary the name for the fun of it
    transform_v2 = ln.Transform(
        description="Introduction v2", revises=transform_v1, version="2"
    ).save()
    assert not transform_v1.is_latest
    assert transform_v2.is_latest
    assert transform_v2.uid.endswith("0001")
    # consciously *not* pass the latest version to revises but the previous
    # it automatically retrieves the latest version
    transform_v3 = ln.Transform(description="Introduction", revises=transform_v1).save()
    assert transform_v3.uid.endswith("0002")
    assert not ln.Transform.objects.get(
        description="Introduction v2", version="2"
    ).is_latest
    assert transform_v3.is_latest
    transform_v4 = ln.Transform(description="Introduction").save()
    assert transform_v4.is_latest
    # add another transform with the same name that's not part of this family
    # but will also be a hit for the query
    assert len(ln.Transform.filter(description="Introduction").all()) == 3
    assert len(ln.Transform.filter(description="Introduction").latest_version()) == 2
    transform_v4.delete()
    assert ln.Transform.get(description="Introduction") == transform_v3
    assert (
        ln.Transform.filter(description="Introduction").latest_version().one()
        == transform_v3
    )

    # test get
    assert ln.Transform.get(transform_v3.uid) == transform_v3
    assert ln.Transform.get(transform_v3.id) == transform_v3
    assert ln.Transform.get(transform_v3.uid[:4]) == transform_v3

    # test delete
    transform_v3.delete()
    assert transform_v2.is_latest

    # test empty QuerySet
    assert (
        ln.Transform.filter(description="IntroductionNotExists")
        .latest_version()
        .one_or_none()
        is None
    )
