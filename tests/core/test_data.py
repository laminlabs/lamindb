import lamindb as ln
import pytest
from lamindb.core._data import add_transform_to_kwargs


def test_add_transform_to_kwargs():
    kwargs = {}
    transform = ln.Transform(name="hello")
    transform.save()
    run = ln.Run(transform)
    add_transform_to_kwargs(kwargs, run)
    assert kwargs["transform"] == transform


def test_rename():
    import pandas as pd
    from lamindb.core.exceptions import RecordNameChangeIntegrityError

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

    curator = ln.Curator.from_df(
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

    # rename label
    ulabel = ln.ULabel.get(name="label-to-rename")
    with pytest.raises(RecordNameChangeIntegrityError):
        ulabel.name = "label-renamed"
        ulabel.save()

    artifact.labels.make_external(ulabel)
    assert not artifact.ulabels.through.objects.filter(
        feature__name="feature_to_rename", ulabel__name="label-to-rename"
    ).exists()
    ulabel.name = "label-renamed"
    ulabel.save()

    # rename feature
    feature = ln.Feature.get(name="feature_to_rename")
    with pytest.raises(RecordNameChangeIntegrityError):
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

    # rename the other feature, automatically deletes no-member featureset
    feature2 = ln.Feature.get(name="feature_to_rename2")
    artifact.features.make_external(feature2)
    assert artifact.feature_sets.count() == 0

    # clean up
    artifact.delete(permanent=True)
    ln.FeatureSet.filter().delete()
    ln.ULabel.filter().delete()
    ln.Feature.filter().delete()
