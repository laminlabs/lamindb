import datetime
import os

import lamindb as ln
import pandas as pd
import pytest


def test_rename_feature(ccaplog):
    df = pd.DataFrame({"old_name": [1, 2]})
    ln.Feature(name="old_name", dtype=int).save()
    artifact = ln.Artifact.from_dataframe(
        df, key="test.parquet", schema="valid_features"
    ).save()
    feature = ln.Feature.get(name="old_name")

    # First rename
    feature.name = "new_name"
    feature.save()
    now1 = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0)
    assert (
        "By renaming feature from 'old_name' to 'new_name' 1 artifact no longer matches the feature name in storage:"
        in ccaplog.text
    )
    if os.getenv("LAMINDB_TEST_DB_VENDOR") != "sqlite":
        feature.refresh_from_db()
        assert feature.synonyms == "old_name"
        assert feature._aux["renamed"] == {
            now1.isoformat().replace("+00:00", "Z"): "old_name"
        }

    # Second rename
    feature.name = "newer_name"
    feature.save()
    now2 = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0)
    assert (
        "By renaming feature from 'new_name' to 'newer_name' 1 artifact no longer matches the feature name in storage:"
        in ccaplog.text
    )
    if os.getenv("LAMINDB_TEST_DB_VENDOR") != "sqlite":
        feature.refresh_from_db()
        assert feature.synonyms == "old_name|new_name"
        assert feature._aux["renamed"] == {
            now1.isoformat().replace("+00:00", "Z"): "old_name",
            now2.isoformat().replace("+00:00", "Z"): "new_name",
        }

    schema = artifact.feature_sets.first()
    artifact.delete(permanent=True)
    schema.delete(permanent=True)
    feature.delete(permanent=True)


@pytest.mark.skipif(
    os.getenv("LAMINDB_TEST_DB_VENDOR") == "sqlite", reason="Postgres-only"
)
def test_rename_label():
    import pandas as pd
    from lamindb.errors import SQLRecordNameChangeIntegrityError

    df = pd.DataFrame(
        {
            "feature_to_rename": [
                "label-to-rename",
                "label-to-rename",
                "label-not-to-rename",
            ],
            "feature_to_rename2": [
                "label-not-to-rename",
                "label-not-to-rename",
                "label-not-to-rename",
            ],
        }
    )

    curator = ln.Curator.from_dataframe(
        df,
        categoricals={
            "feature_to_rename": ln.ULabel.name,
            "feature_to_rename2": ln.ULabel.name,
        },
    )
    curator.add_new_from("feature_to_rename")
    curator.add_new_from("feature_to_rename2")
    artifact = curator.save_artifact(description="test-rename")
    assert artifact.ulabels.through.objects.filter(
        feature__name="feature_to_rename", ulabel__name="label-to-rename"
    ).exists()
    assert ln.Artifact.filter(feature_sets__features__name="feature_to_rename").exists()

    # rename feature
    feature = ln.Feature.get(name="feature_to_rename")
    with pytest.raises(SQLRecordNameChangeIntegrityError):
        feature.name = "feature_renamed"
        feature.save()

    artifact.features.make_external(feature)
    assert not ln.Artifact.filter(
        feature_sets__features__name="feature_to_rename"
    ).exists()
    assert ln.Artifact.filter(
        feature_sets__features__name="feature_to_rename2"
    ).exists()
    feature.name = "feature_renamed"
    feature.save()

    # rename the other feature, automatically deletes no-member schema
    feature2 = ln.Feature.get(name="feature_to_rename2")
    artifact.features.make_external(feature2)
    assert artifact.feature_sets.count() == 0

    # clean up
    artifact.delete(permanent=True)
    ln.Schema.filter().delete(permanent=True)
    ln.ULabel.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
