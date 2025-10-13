import lamindb as ln
import pandas as pd
import pytest
from lamindb.models._is_versioned import bump_version, set_version


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


def test_add_to_version_family(df1, df2):
    artifact1 = ln.Artifact.from_dataframe(df1, description="test1")
    artifact1.save()
    artifact2 = ln.Artifact.from_dataframe(df2, description="test2")
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


def test_transform_versioning_based_on_key():
    transform1 = ln.Transform(
        key="test-pipeline",
        version="1.0",
        source_code="1",
        type="pipeline",
    ).save()
    assert transform1.is_latest

    with pytest.raises(ValueError) as e:
        transform2 = ln.Transform(
            key="test-pipeline",
            version="1.0",
            source_code="2",
            type="pipeline",
        ).save()
    assert (
        e.exconly()
        == "ValueError: Please change the version tag or leave it `None`, '1.0' is already taken"
    )

    transform2 = ln.Transform(
        key="test-pipeline",
        # do not pass the version tag, which corresponds to: version=None
        source_code="2",
        type="pipeline",
    ).save()

    assert transform2.version is None
    assert transform2.is_latest
    assert transform2.hash != transform1.hash
    assert not ln.Transform.get(key="test-pipeline", version="1.0").is_latest

    transform3 = ln.Transform(
        key="test-pipeline",
        version="abcd",  # mimic commit hash
        source_code="3",
        type="pipeline",
    ).save()

    assert transform3.version == "abcd"
    assert transform3.is_latest
    assert transform3.hash != transform2.hash
    assert not ln.Transform.get(key="test-pipeline", source_code="2").is_latest


def test_transform_versioning_based_on_revises():
    # build one version family
    transform_v1 = ln.Transform(description="Introduction").save()
    assert transform_v1.is_latest
    assert transform_v1.version is None

    # pass the latest version
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

    # add another transform with the same description that's not part of this family
    # but will also be a hit for the query
    assert len(ln.Transform.filter(description="Introduction")) == 3
    assert len(ln.Transform.filter(description="Introduction").latest_version()) == 2
    transform_v4.delete(permanent=True)
    assert ln.Transform.get(description="Introduction") == transform_v3
    assert (
        ln.Transform.filter(description="Introduction").latest_version().one()
        == transform_v3
    )

    # test get
    assert ln.Transform.get(transform_v3.uid) == transform_v3
    assert ln.Transform.get(transform_v3.id) == transform_v3
    assert ln.Transform.get(transform_v3.uid[:-4]) == transform_v3

    # test empty QuerySet
    assert (
        ln.Transform.filter(description="IntroductionNotExists")
        .latest_version()
        .one_or_none()
        is None
    )

    # test soft delete
    transform_v3.delete()
    assert transform_v2.is_latest

    # test hard delete
    transform_v2.delete(permanent=True)
    assert (
        transform_v1_retrieved := ln.Transform.get(transform_v3.uid[:-4])
    ) == transform_v1
    assert transform_v1_retrieved.is_latest

    # test soft delete on the last existing version does not change is_latest
    transform_v1_retrieved.delete()
    assert (
        transform_v1_retrieved := ln.Transform.get(transform_v1.uid)
    ) == transform_v1
    assert transform_v1_retrieved.is_latest

    # fully delete
    transform_v1.delete(permanent=True)

    # last object that exists is in the trash
    assert ln.Transform.get(transform_v3.uid[:-4]) == transform_v3
    assert transform_v3.branch_id == -1
    transform_v3.delete(permanent=True)


def test_path_rename():
    # this is related to renames inside _add_to_version_family
    with open("test_new_path.txt", "w") as f:
        f.write("test_new_path")
    old_path = ln.UPath("s3://lamindata/.lamindb/test_new_path.txt")
    old_path.upload_from("./test_new_path.txt")
    assert old_path.exists()
    new_path = old_path.rename(old_path.with_name("test_new_path2.txt"))
    assert new_path.exists()
    assert new_path.as_posix() == "s3://lamindata/.lamindb/test_new_path2.txt"
    assert not old_path.exists()
    new_path.unlink()
    ln.UPath("./test_new_path.txt").unlink()
