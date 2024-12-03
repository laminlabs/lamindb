# .latest_version is tested in test_versioning.py

import bionty as bt
import lamindb as ln
import pytest
from lamindb._query_set import DoesNotExist
from lnschema_core.users import current_user_id


# please also see the test_curate_df.py tests
def test_df():
    project_label = ln.ULabel(name="project").save()
    project_names = [f"Project {i}" for i in range(3)]
    labels = ln.ULabel.from_values(project_names, create=True).save()
    project_label.children.add(*labels)
    df = ln.ULabel.df(include="parents__name")
    assert df.columns[3] == "parents__name"
    assert df["parents__name"].iloc[0] == {project_label.name}
    df = ln.ULabel.df(include=["parents__name", "parents__created_by_id"])
    assert df.columns[4] == "parents__created_by_id"
    assert df["parents__name"].iloc[0] == {project_label.name}
    assert set(df["parents__created_by_id"].iloc[0]) == {current_user_id()}

    # for other models
    feature_names = [f"Feature {i}" for i in range(3)]
    features = [ln.Feature(name=name, dtype=int) for name in feature_names]
    ln.save(features)
    feature_set = ln.FeatureSet(features, name="my feature_set").save()
    feature_set.features.set(features)

    df = ln.FeatureSet.filter(name="my feature_set").df(include="features__name")
    assert df.columns[3] == "features__name"
    # order is not conserved
    assert set(df["features__name"].iloc[0]) == set(feature_names)
    # pass a list
    df = ln.FeatureSet.filter(name="my feature_set").df(
        include=["features__name", "features__created_by_id"]
    )
    assert df.columns[4] == "features__created_by_id"
    assert set(df["features__name"].iloc[0]) == set(feature_names)
    assert set(df["features__created_by_id"].iloc[0]) == {current_user_id()}

    # inner join parents on features
    df = ln.FeatureSet.filter().df(
        include=["features__name", "features__created_by_id"]
    )
    assert set(df["features__name"].iloc[0]) == set(feature_names)
    assert set(df["features__created_by_id"].iloc[0]) == {current_user_id()}

    # raise error for non many-to-many
    df = ln.ULabel.filter(name="Project 0").df(include="created_by__name")
    assert df["created_by__name"].iloc[0] == "Test User1"

    # do not return fields with no data in the registry
    # does not make sense in Alex's opinion
    # too much magic; got removed in https://github.com/laminlabs/lamindb/pull/2238
    # df = (
    #     ln.Artifact.using("laminlabs/cellxgene")
    #     .filter(suffix=".h5ad")
    #     .df(include=["tissues__name", "pathways__name"])
    # )
    # assert "tissues__name" in df.columns
    # assert "pathways__name" not in df.columns
    # assert df.shape[0] > 0

    # clean up
    project_label.delete()
    for label in labels:
        label.delete()

    feature_set.delete()
    for feature in features:
        feature.delete()

    # call it from a non-select-derived queryset
    qs = ln.User.objects.all()
    assert qs.df().iloc[0]["handle"] == "testuser1"


def test_one_first():
    qs = ln.User.objects.all()
    assert qs.one().handle == "testuser1"
    assert qs.first().handle == "testuser1"
    assert qs.one_or_none().handle == "testuser1"

    qs = ln.User.filter(handle="test")
    with pytest.raises(DoesNotExist):
        qs.one()
    qs = bt.Source.filter().all()
    # should be MultipleResultsFound, but internal to Django
    with pytest.raises(Exception):  # noqa: B017
        qs.one()
    with pytest.raises(Exception):  # noqa: B017
        qs.one_or_none()


def test_search():
    label_names = [f"ULabel {i}" for i in range(3)]
    labels = [ln.ULabel(name=name) for name in label_names]
    ln.save(labels)
    qs = ln.ULabel.filter(name__startswith="ULabel").all()
    assert qs.search("ULabel 1")[0].name == "ULabel 1"
    assert qs.search("ULabel 1", field=ln.ULabel.name)[0].name == "ULabel 1"
    for label in labels:
        label.delete()


def test_lookup():
    qs = ln.User.filter(handle="testuser1").all()
    # pass str to field
    lookup = qs.lookup(field="handle")
    assert lookup.testuser1.handle == "testuser1"
    # pass StrField to field
    lookup = qs.lookup(field=ln.User.handle)
    assert lookup.testuser1.handle == "testuser1"
    # manager, default field
    qsm = ln.User.filter(handle="testuser1")
    lookup = qsm.lookup()
    assert lookup.testuser1.handle == "testuser1"


def test_inspect():
    qs = ln.User.filter(handle="testuser1").all()
    assert qs.inspect(["user1", "user2"], "name")["validated"] == []
    assert ln.User.inspect(["user1", "user2"], "name")["validated"] == []
    assert ln.User.inspect(["user1", "user2"], ln.User.name)["validated"] == []
    assert ln.User.inspect("user1", "name")["validated"] == []


def test_validate():
    qs = ln.User.filter(handle="testuser1").all()
    assert qs.validate(["testuser1", "Test User1"], "handle").tolist() == [True, False]
    assert ln.User.validate(["testuser1", "Test User1"], "handle").tolist() == [
        True,
        False,
    ]
    assert ln.User.validate(["testuser1", "Test User1"], ln.User.handle).tolist() == [
        True,
        False,
    ]
    # returns True
    assert ln.User.validate("testuser1", ln.User.handle)


def test_map_synonyms():
    qs = ln.User.filter(handle="testuser1").all()
    assert qs.standardize(["user1", "user2"]) == ["user1", "user2"]
