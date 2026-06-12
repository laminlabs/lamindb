import lamindb as ln
import pandas as pd
import pytest
from lamindb.errors import IntegrityError
from lamindb.models._is_versioned import (
    _adjust_is_latest_when_deleting_is_versioned,
    bump_version,
    max_version_uid_in_family,
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


def test_transform_versioning_across_branches_preserves_main_latest():
    main_branch = ln.Branch.get(name="main")
    ln.setup.switch(main_branch.name)
    branch = ln.Branch(name="test_versioning_branch_latest").save()
    transform_v1 = ln.Transform(
        key="test-branch-aware-is-latest",
        source_code="main-v1",
        kind="pipeline",
    ).save()
    try:
        ln.setup.switch(branch.name)
        transform_v2 = ln.Transform(
            key="test-branch-aware-is-latest",
            revises=transform_v1,
            source_code="feature-v2",
            kind="pipeline",
        ).save()
        transform_v1.refresh_from_db()
        assert transform_v1.is_latest
        assert transform_v2.is_latest

        # Passing an older revises still increments from the family max uid.
        transform_v3 = ln.Transform(
            key="test-branch-aware-is-latest",
            revises=transform_v1,
            source_code="feature-v3",
            kind="pipeline",
        ).save()
        transform_v2.refresh_from_db()
        transform_v1.refresh_from_db()
        assert transform_v3.uid.endswith("0002")
        assert not transform_v2.is_latest
        assert transform_v3.is_latest
        assert transform_v1.is_latest
    finally:
        ln.setup.switch(main_branch.name)
        for uid in (transform_v1.uid[:-4],):
            for record in ln.Transform.objects.filter(uid__startswith=uid):
                record.delete(permanent=True)
        branch.delete(permanent=True)


def test_artifact_versioning_across_branches_preserves_main_latest():
    main_branch = ln.Branch.get(name="main")
    ln.setup.switch(main_branch.name)
    branch = ln.Branch(name="test_artifact_versioning_branch_latest").save()
    artifact_v1 = ln.Artifact.from_dataframe(
        pd.DataFrame({"branch_feat": [10, 11]}),
        key="test-artifact-branch-aware-is-latest.parquet",
        description="main-v1",
    ).save()
    # sanity check: the artifact is created on the branch steered by `switch`.
    assert artifact_v1.branch_id == main_branch.id
    try:
        ln.setup.switch(branch.name)
        artifact_v2 = ln.Artifact.from_dataframe(
            pd.DataFrame({"branch_feat": [12, 13]}),
            revises=artifact_v1,
            description="feature-v2",
        ).save()
        assert artifact_v2.branch_id == branch.id
        artifact_v1.refresh_from_db()
        # main's head stays latest; the contribution branch gets its own head.
        assert artifact_v1.is_latest
        assert artifact_v2.is_latest
        assert artifact_v2.uid.endswith("0001")

        # Passing an older `revises` still increments from the family max uid and only
        # demotes the head on the *creation* branch, leaving main's head intact.
        artifact_v3 = ln.Artifact.from_dataframe(
            pd.DataFrame({"branch_feat": [14, 15]}),
            revises=artifact_v1,
            description="feature-v3",
        ).save()
        assert artifact_v3.branch_id == branch.id
        artifact_v2.refresh_from_db()
        artifact_v1.refresh_from_db()
        assert artifact_v3.uid.endswith("0002")
        assert artifact_v3.stem_uid == artifact_v1.stem_uid
        assert not artifact_v2.is_latest
        assert artifact_v3.is_latest
        assert artifact_v1.is_latest
    finally:
        ln.setup.switch(main_branch.name)
        for record in ln.Artifact.objects.filter(uid__startswith=artifact_v1.uid[:-4]):
            record.delete(permanent=True)
        branch.delete(permanent=True)


def test_stale_revises_raises_integrity_error():
    transform_v1 = ln.Transform(
        key="stale-revises-validation-error",
        source_code="v1",
        kind="pipeline",
    ).save()
    transform_v2 = ln.Transform(
        key="stale-revises-validation-error",
        revises=transform_v1,
        source_code="v2",
        kind="pipeline",
    ).save()
    try:
        transform_pending = ln.Transform(
            key="stale-revises-validation-error",
            revises=transform_v2,
            source_code="v3",
            kind="pipeline",
        )
        assert not transform_pending._refresh_revises_if_stale
        # Simulate a concurrent writer demoting transform_v2 after init but before save.
        ln.Transform(
            key="stale-revises-validation-error",
            revises=transform_v2,
            source_code="v4",
            kind="pipeline",
        ).save()

        with pytest.raises(IntegrityError) as error:
            transform_pending.save()
        message = str(error.value)
        assert "Cannot revise a non-latest record" in message
        assert f"revises=Transform(uid={transform_v2.uid}" in message
        assert "key=stale-revises-validation-error" in message
        assert "new=Transform(uid=" in message
    finally:
        for record in ln.Transform.filter(uid__startswith=transform_v1.stem_uid):
            record.delete(permanent=True)


def test_inferred_revises_refreshes_and_requeries_latest():
    transform_v1 = ln.Transform(
        key="stale-inferred-revises-requery",
        source_code="v1",
        kind="pipeline",
    ).save()
    transform_v2 = ln.Transform(
        key="stale-inferred-revises-requery",
        revises=transform_v1,
        source_code="v2",
        kind="pipeline",
    ).save()
    try:
        transform_pending = ln.Transform(
            key="stale-inferred-revises-requery",
            source_code="v3",
            kind="pipeline",
        )
        assert transform_pending._revises is not None
        assert transform_pending._revises.uid == transform_v2.uid
        assert transform_pending._refresh_revises_if_stale

        # Simulate stale latest flags without creating a new version.
        ln.Transform.filter(id=transform_v2.id).update(is_latest=False)
        ln.Transform.filter(id=transform_v1.id).update(is_latest=True)

        # If refresh/requery is disabled, save should fail on stale revises.
        transform_pending._refresh_revises_if_stale = False
        with pytest.raises(IntegrityError) as error:
            transform_pending.save()
        assert "Cannot revise a non-latest record" in str(error.value)

        # Re-enable behavior and verify save now succeeds.
        transform_pending._refresh_revises_if_stale = True
        transform_pending.save()
        assert transform_pending.is_latest
        transform_pending.refresh_from_db()
        assert transform_pending.is_latest
    finally:
        for record in ln.Transform.filter(uid__startswith=transform_v1.stem_uid):
            record.delete(permanent=True)


def test_inferred_revises_handles_concurrent_new_version():
    transform_v1 = ln.Transform(
        key="stale-inferred-revises-concurrent",
        source_code="v1",
        kind="pipeline",
    ).save()
    transform_v2 = ln.Transform(
        key="stale-inferred-revises-concurrent",
        revises=transform_v1,
        source_code="v2",
        kind="pipeline",
    ).save()
    try:
        # inferred revises picks the current latest (v2); uid is provisionally ...0002
        transform_pending = ln.Transform(
            key="stale-inferred-revises-concurrent",
            source_code="v3",
            kind="pipeline",
        )
        assert transform_pending._refresh_revises_if_stale
        assert transform_pending._revises.uid == transform_v2.uid
        assert transform_pending.uid.endswith("0002")

        # Simulate a concurrent writer creating a genuinely new version (new uid)
        # after `transform_pending` was initialized but before it is saved. This
        # advances the family max uid to ...0002, which `transform_pending` would
        # collide with if its uid were not recomputed at save time.
        transform_concurrent = ln.Transform(
            key="stale-inferred-revises-concurrent",
            revises=transform_v2,
            source_code="v-concurrent",
            kind="pipeline",
        ).save()
        assert transform_concurrent.uid.endswith("0002")

        # Saving must create the next version (...0003) from the new family head and
        # keep this record's own content, not collide with / silently return the
        # concurrently-created record.
        transform_pending.save()
        assert transform_pending.uid.endswith("0003")
        assert transform_pending.source_code == "v3"
        assert transform_pending.is_latest
        transform_pending.refresh_from_db()
        assert transform_pending.is_latest
        transform_concurrent.refresh_from_db()
        assert not transform_concurrent.is_latest
    finally:
        for record in ln.Transform.filter(uid__startswith=transform_v1.stem_uid):
            record.delete(permanent=True)


def test_version_chain_crosses_base62_case_boundary():
    # base62 increments digits -> uppercase -> lowercase, so the version-suffix chain
    # crosses ...000Z -> ...000a -> ...000b. This must hold regardless of the database
    # collation: locale-aware collations (e.g. Postgres `en_US.UTF-8`) sort the letter
    # `Z` after `a`, which previously made the chain stall at the Z -> a boundary
    # (`...000Z` kept being selected as the latest, regenerating `...000a` forever).
    transform_v1 = ln.Transform(
        key="version-chain-z-to-a",
        source_code="v1",
        kind="pipeline",
    ).save()
    stem_uid = transform_v1.stem_uid
    try:
        # fast-forward the version suffix to ...000Z without creating 35 versions
        ln.Transform.filter(id=transform_v1.id).update(uid=stem_uid + "000Z")
        transform_v1 = ln.Transform.get(id=transform_v1.id)
        assert transform_v1.uid.endswith("000Z")
        assert transform_v1.is_latest

        # crossing the Z -> a boundary
        transform_v2 = ln.Transform(
            key="version-chain-z-to-a",
            revises=transform_v1,
            source_code="v2",
            kind="pipeline",
        ).save()
        assert transform_v2.uid.endswith("000a")

        # the family now holds both ...000Z and ...000a; the max must be ...000a
        # (base62), not ...000Z (locale collation) -- assert the selector directly.
        assert max_version_uid_in_family(transform_v2).endswith("000a")

        # so the next version is ...000b rather than colliding back on ...000a
        transform_v3 = ln.Transform(
            key="version-chain-z-to-a",
            revises=transform_v2,
            source_code="v3",
            kind="pipeline",
        ).save()
        assert transform_v3.uid.endswith("000b")
        assert transform_v3.is_latest
        transform_v1.refresh_from_db()
        transform_v2.refresh_from_db()
        assert not transform_v1.is_latest
        assert not transform_v2.is_latest
    finally:
        for record in ln.Transform.filter(uid__startswith=stem_uid):
            record.delete(permanent=True)


def test_inferred_revises_prefers_target_branch_head():
    main_branch = ln.Branch.get(name="main")
    ln.setup.switch(main_branch.name)
    branch = ln.Branch(name="test_inferred_revises_branch").save()
    transform_main_v1 = ln.Transform(
        key="inferred-revises-branch-aware",
        source_code="main-v1",
        kind="pipeline",
    ).save()
    # sanity check: `switch` must actually steer the creation branch, otherwise the
    # rest of the test (which relies on heads living on distinct branches) is moot.
    assert transform_main_v1.branch_id == main_branch.id
    try:
        # create a *more recently* created head on another branch for the same key
        ln.setup.switch(branch.name)
        transform_branch = ln.Transform(
            key="inferred-revises-branch-aware",
            source_code="branch-v1",
            kind="pipeline",
        ).save()
        assert transform_branch.branch_id == branch.id
        assert transform_branch.is_latest
        transform_main_v1.refresh_from_db()
        # main head is preserved (cross-branch revision does not demote it)
        assert transform_main_v1.is_latest

        # back on main, infer revises (no explicit `revises`). Even though the other
        # branch's head was created more recently, the inferred revises must be the
        # head on the *target* (main) branch.
        ln.setup.switch(main_branch.name)
        transform_main_v2 = ln.Transform(
            key="inferred-revises-branch-aware",
            source_code="main-v2",
            kind="pipeline",
        )
        # the record-to-be must target main, so the inferred revises is scoped to main
        assert transform_main_v2.branch_id == main_branch.id
        assert transform_main_v2._revises is not None
        assert transform_main_v2._revises.uid == transform_main_v1.uid
        transform_main_v2.save()

        transform_main_v1.refresh_from_db()
        transform_branch.refresh_from_db()
        assert transform_main_v2.branch_id == main_branch.id
        assert transform_main_v2.is_latest
        # the main head was demoted, the other branch's head is untouched
        assert not transform_main_v1.is_latest
        assert transform_branch.is_latest
        # exactly one latest on main for this family (no double head)
        main_latest = ln.Transform.objects.filter(
            uid__startswith=transform_main_v2.stem_uid,
            branch_id=main_branch.id,
            is_latest=True,
        )
        assert main_latest.count() == 1
    finally:
        ln.setup.switch(main_branch.name)
        for record in ln.Transform.objects.filter(
            uid__startswith=transform_main_v1.stem_uid
        ):
            record.delete(permanent=True)
        branch.delete(permanent=True)


def test_inferred_revises_prefers_target_branch_head_artifact():
    main_branch = ln.Branch.get(name="main")
    ln.setup.switch(main_branch.name)
    branch = ln.Branch(name="test_inferred_revises_branch_artifact").save()
    key = "inferred-revises-branch-aware-artifact.parquet"
    artifact_main_v1 = ln.Artifact.from_dataframe(
        pd.DataFrame({"inferred_feat": [20, 21]}),
        key=key,
        description="main-v1",
    ).save()
    # sanity check: `switch` must actually steer the creation branch, otherwise the
    # rest of the test (which relies on heads living on distinct branches) is moot.
    assert artifact_main_v1.branch_id == main_branch.id
    try:
        # create a *more recently* created head on another branch for the same key;
        # with no explicit `revises` it infers the family head (here only main's) and
        # joins the family on the contribution branch, leaving main's head intact.
        ln.setup.switch(branch.name)
        artifact_branch = ln.Artifact.from_dataframe(
            pd.DataFrame({"inferred_feat": [22, 23]}),
            key=key,
            description="branch-v1",
        ).save()
        assert artifact_branch.branch_id == branch.id
        assert artifact_branch.is_latest
        assert artifact_branch.stem_uid == artifact_main_v1.stem_uid
        artifact_main_v1.refresh_from_db()
        # main head is preserved (cross-branch revision does not demote it)
        assert artifact_main_v1.is_latest

        # back on main, infer revises (no explicit `revises`). Even though the other
        # branch's head was created more recently, the inferred revises must be the
        # head on the *target* (main) branch.
        ln.setup.switch(main_branch.name)
        artifact_main_v2 = ln.Artifact.from_dataframe(
            pd.DataFrame({"inferred_feat": [24, 25]}),
            key=key,
            description="main-v2",
        )
        # the record-to-be must target main, so the inferred revises is scoped to main
        assert artifact_main_v2.branch_id == main_branch.id
        assert artifact_main_v2._revises is not None
        assert artifact_main_v2._revises.uid == artifact_main_v1.uid
        artifact_main_v2.save()

        artifact_main_v1.refresh_from_db()
        artifact_branch.refresh_from_db()
        assert artifact_main_v2.branch_id == main_branch.id
        assert artifact_main_v2.is_latest
        # the main head was demoted, the other branch's head is untouched
        assert not artifact_main_v1.is_latest
        assert artifact_branch.is_latest
        # exactly one latest on main for this family (no double head)
        main_latest = ln.Artifact.objects.filter(
            uid__startswith=artifact_main_v2.stem_uid,
            branch_id=main_branch.id,
            is_latest=True,
        )
        assert main_latest.count() == 1
    finally:
        ln.setup.switch(main_branch.name)
        for record in ln.Artifact.objects.filter(
            uid__startswith=artifact_main_v1.stem_uid
        ):
            record.delete(permanent=True)
        branch.delete(permanent=True)


def test_soft_delete_promotes_within_branch():
    main_branch = ln.Branch.get(name="main")
    ln.setup.switch(main_branch.name)
    branch = ln.Branch(name="test_soft_delete_within_branch").save()
    key = "soft-delete-branch-aware.parquet"
    # a version family with two versions on the contribution branch
    ln.setup.switch(branch.name)
    artifact_b1 = ln.Artifact.from_dataframe(
        pd.DataFrame({"sd_feat": [1, 2]}), key=key, description="b1"
    ).save()
    artifact_b2 = ln.Artifact.from_dataframe(
        pd.DataFrame({"sd_feat": [3, 4]}), revises=artifact_b1, description="b2"
    ).save()
    artifact_b1.refresh_from_db()
    assert not artifact_b1.is_latest
    assert artifact_b2.is_latest
    try:
        # a *more recently* created head on main for the same family
        ln.setup.switch(main_branch.name)
        artifact_main = ln.Artifact.from_dataframe(
            pd.DataFrame({"sd_feat": [5, 6]}), revises=artifact_b1, description="main"
        ).save()
        assert artifact_main.branch_id == main_branch.id
        assert artifact_main.is_latest
        artifact_b2.refresh_from_db()
        # cross-branch revise leaves the contribution branch's head intact
        assert artifact_b2.is_latest

        # soft-delete the branch head: the next version *on the branch* must be
        # promoted, not main's (more recently created) head.
        ln.setup.switch(branch.name)
        artifact_b2.delete()
        artifact_b1.refresh_from_db()
        artifact_main.refresh_from_db()
        assert artifact_b1.is_latest  # promoted within the branch
        assert artifact_b1.branch_id == branch.id
        assert artifact_main.is_latest  # main head untouched
        # exactly one latest per branch for this family (no headless branch, no double head)
        for branch_id in (branch.id, main_branch.id):
            assert (
                ln.Artifact.objects.filter(
                    uid__startswith=artifact_b1.stem_uid,
                    branch_id=branch_id,
                    is_latest=True,
                ).count()
                == 1
            )
    finally:
        ln.setup.switch(main_branch.name)
        for record in ln.Artifact.objects.filter(uid__startswith=artifact_b1.stem_uid):
            record.delete(permanent=True)
        branch.delete(permanent=True)


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
    """Direct unit test for _adjust_is_latest_when_deleting_is_versioned (covers multiple promoted)."""
    # Build two version families, each with v1 (older) and v2 (latest)
    v1a = ln.Transform(key="Adjust latest family A").save()
    v2a = ln.Transform(revises=v1a, key="Adjust latest family A").save()
    v1b = ln.Transform(key="Adjust latest family B").save()
    v2b = ln.Transform(revises=v1b, key="Adjust latest family B").save()
    assert v2a.is_latest and v2b.is_latest
    assert not v1a.is_latest and not v1b.is_latest

    # Delete both latest → two promoted (covers "new latest ... versions: [...]" branch)
    promoted = _adjust_is_latest_when_deleting_is_versioned([v2a, v2b])
    assert len(promoted) == 2
    assert set(promoted) == {v1a.pk, v1b.pk}

    v1a.refresh_from_db()
    v1b.refresh_from_db()
    assert v1a.is_latest and v1b.is_latest

    # Edge case: empty list returns []
    assert _adjust_is_latest_when_deleting_is_versioned([]) == []

    # Clean up
    v2a.delete(permanent=True)
    v2b.delete(permanent=True)
    v1a.delete(permanent=True)
    v1b.delete(permanent=True)
