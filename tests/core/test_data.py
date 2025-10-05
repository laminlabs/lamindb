import lamindb as ln
import pytest


def test_rename():
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
            "feature_to_rename": ln.Record.name,
            "feature_to_rename2": ln.Record.name,
        },
    )
    curator.add_new_from("feature_to_rename")
    curator.add_new_from("feature_to_rename2")
    artifact = curator.save_artifact(description="test-rename")
    assert artifact.records.through.objects.filter(
        feature__name="feature_to_rename", record__name="label-to-rename"
    ).exists()
    assert ln.Artifact.filter(feature_sets__features__name="feature_to_rename").exists()

    # rename label
    record = ln.Record.get(name="label-to-rename")
    with pytest.raises(SQLRecordNameChangeIntegrityError):
        record.name = "label-renamed"
        record.save()

    artifact.labels.make_external(record)
    assert not artifact.records.through.objects.filter(
        feature__name="feature_to_rename", record__name="label-to-rename"
    ).exists()
    record.name = "label-renamed"
    record.save()

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
    ln.Record.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
