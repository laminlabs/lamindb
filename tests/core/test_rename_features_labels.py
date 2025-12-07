import datetime
import os

import lamindb as ln
import pandas as pd


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
        "by renaming feature from 'old_name' to 'new_name' 1 artifact no longer matches the feature name in storage:"
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
        "by renaming feature from 'new_name' to 'newer_name' 1 artifact no longer matches the feature name in storage:"
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


def test_rename_label(ccaplog):
    df = pd.DataFrame(
        {
            "feature1": pd.Categorical(["label1", "label2"]),
            "feature2": pd.Categorical(["label2", "label2"]),
        }
    )

    ulabel1 = ln.ULabel(name="label1").save()
    ulabel2 = ln.ULabel(name="label2").save()
    feature1 = ln.Feature(name="feature1", dtype=ln.ULabel).save()
    feature2 = ln.Feature(name="feature2", dtype=ln.ULabel).save()
    artifact = ln.Artifact.from_dataframe(
        df, key="test.parquet", schema="valid_features"
    ).save()

    ulabel = ln.ULabel.get(name="label1")
    ulabel.name = "label-renamed"
    ulabel.save()

    assert (
        "by renaming label from 'label1' to 'label-renamed' 1 artifact no longer matches the label name in storage:"
        in ccaplog.text
    )

    schema = artifact.feature_sets.first()
    artifact.delete(permanent=True)
    schema.delete(permanent=True)
    feature1.delete(permanent=True)
    feature2.delete(permanent=True)
    ulabel1.delete(permanent=True)
    ulabel2.delete(permanent=True)
