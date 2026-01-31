import lamindb as ln
import pandas as pd
import pytest
from lamindb.models._is_versioned import (
    _adjust_is_latest_when_deleting_is_versioned,
    bump_version,
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


def test_add_to_version_family(df1, df2):
    artifact1 = ln.Artifact.from_dataframe(df1, description="test1").save()
    artifact2 = ln.Artifact.from_dataframe(df2, description="test2").save()
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
        kind="pipeline",
    ).save()
    assert transform1.is_latest
    assert transform1.version_tag == "1.0"
    assert transform1.version == "1.0"

    with pytest.raises(ValueError) as e:
        transform2 = ln.Transform(
            key="test-pipeline",
            version="1.0",
            source_code="2",
            kind="pipeline",
        ).save()
    assert (
        e.exconly()
        == "ValueError: Please change the version tag or leave it `None`, '1.0' is already taken"
    )

    transform2 = ln.Transform(
        key="test-pipeline",
        # do not pass the version tag, which corresponds to: version=None
        source_code="2",
        kind="pipeline",
    ).save()

    assert transform2.version_tag is None
    assert transform2.version == transform2.uid[-4:]  # version falls back to uid suffix
    assert transform2.is_latest
    assert transform2.hash != transform1.hash
    assert not ln.Transform.get(key="test-pipeline", version="1.0").is_latest

    transform3 = ln.Transform(
        key="test-pipeline",
        version="abcd",  # mimic commit hash
        source_code="3",
        kind="pipeline",
    ).save()

    assert transform3.version_tag == "abcd"
    assert transform3.version == "abcd"
    assert transform3.is_latest
    assert transform3.hash != transform2.hash
    assert not ln.Transform.get(key="test-pipeline", source_code="2").is_latest


def test_transform_versioning_based_on_revises():
    # build one version family
    transform_v1 = ln.Transform(key="Introduction").save()
    assert transform_v1.is_latest
    assert transform_v1.version_tag is None

    # pass the latest version
    transform_v2 = ln.Transform(
        key="Introduction v2", revises=transform_v1, version="2"
    ).save()
    assert not transform_v1.is_latest
    assert transform_v2.is_latest
    assert transform_v2.uid.endswith("0001")
    assert transform_v2.version_tag == "2"
    assert transform_v2.version == "2"

    # consciously *not* pass the latest version to revises but the previous
    # it automatically retrieves the latest version
    transform_v3 = ln.Transform(key="Introduction", revises=transform_v1).save()
    assert transform_v3.uid.endswith("0002")
    assert not ln.Transform.get(key="Introduction v2", version="2").is_latest
    assert transform_v3.is_latest
    # no source code code was yet saved, returning existing transform with same key
    transform_v4 = ln.Transform(key="Introduction").save()
    assert transform_v4 == transform_v3

    assert len(ln.Transform.filter(key="Introduction")) == 2
    assert len(ln.Transform.filter(key="Introduction").filter(is_latest=True)) == 1
    assert ln.Transform.get(key="Introduction") == transform_v3
    assert ln.Transform.filter(key="Introduction").get(is_latest=True) == transform_v3

    # test get
    assert ln.Transform.get(transform_v3.uid) == transform_v3
    assert ln.Transform.get(transform_v3.id) == transform_v3
    assert ln.Transform.get(transform_v3.uid[:-4]) == transform_v3

    # test empty QuerySet
    assert (
        ln.Transform.filter(key="IntroductionNotExists")
        .filter(is_latest=True)
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


def test_version_backward_compatibility():
    """Test that queries using version= still work (backward compatibility)."""
    # Create transforms with different versions and source_code to avoid deduplication
    transform1 = ln.Transform(
        key="test-backward-compat",
        version="1.0",
        kind="pipeline",
        source_code="code1",
    ).save()
    transform2 = ln.Transform(
        key="test-backward-compat",
        version="2.0",
        kind="pipeline",
        source_code="code2",
    ).save()

    # Test that we can query using version= (old API)
    found = ln.Transform.get(key="test-backward-compat", version="1.0")
    assert found == transform1
    assert found.version_tag == "1.0"
    assert found.version == "1.0"

    found = ln.Transform.get(key="test-backward-compat", version="2.0")
    assert found == transform2
    assert found.version_tag == "2.0"
    assert found.version == "2.0"

    # Test filter with version=
    results = ln.Transform.filter(key="test-backward-compat", version="1.0")
    assert len(results) == 1
    assert results.first() == transform1

    # Test with Artifact
    artifact1 = ln.Artifact.from_dataframe(
        pd.DataFrame({"col1": [1, 2]}), key="test-artifact.parquet", version="1.0"
    ).save()
    artifact2 = ln.Artifact.from_dataframe(
        pd.DataFrame({"col1": [3, 4]}), key="test-artifact.parquet", version="2.0"
    ).save()

    found_artifact = ln.Artifact.get(key="test-artifact.parquet", version="1.0")
    assert found_artifact == artifact1
    assert found_artifact.version_tag == "1.0"
    assert found_artifact.version == "1.0"

    found_artifact = ln.Artifact.get(key="test-artifact.parquet", version="2.0")
    assert found_artifact == artifact2
    assert found_artifact.version_tag == "2.0"
    assert found_artifact.version == "2.0"

    # Cleanup
    transform1.delete(permanent=True)
    transform2.delete(permanent=True)
    artifact1.delete(permanent=True)
    artifact2.delete(permanent=True)


def test_adjust_is_latest_when_deleting_is_versioned():
    """Direct unit test for _adjust_is_latest_when_deleting_is_versioned."""
    # Build one version family: v1 (older), v2 (latest)
    v1 = ln.Transform(key="Adjust latest unit test").save()
    v2 = ln.Transform(revises=v1, key="Adjust latest unit test").save()
    assert v2.is_latest
    assert not v1.is_latest

    db = getattr(v1._state, "db", None) or "default"
    promoted = _adjust_is_latest_when_deleting_is_versioned(ln.Transform, db, [v2.pk])
    assert len(promoted) == 1
    assert promoted[0].pk == v1.pk

    v1.refresh_from_db()
    assert v1.is_latest

    # Edge case: empty id_list returns []
    assert _adjust_is_latest_when_deleting_is_versioned(ln.Transform, db, []) == []

    # Clean up (v2 first so v1 stays sole latest, then v1)
    v2.delete(permanent=True)
    v1.delete(permanent=True)
