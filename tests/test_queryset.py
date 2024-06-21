# .latest_version is tested in test_versioning.py

import bionty as bt
import lamindb as ln
import pytest
from lamindb._query_set import MultipleResultsFound, NoResultFound
from lnschema_core.users import current_user_id


def test_df():
    # for self-referential models
    project_label = ln.ULabel(name="Project")
    project_label.save()
    project_names = [f"Project {i}" for i in range(3)]
    labels = [ln.ULabel(name=name) for name in project_names]
    ln.save(labels)
    for label in labels:
        label.parents.add(project_label)
    df = ln.ULabel.filter().df(include="parents__name")
    assert df.columns[0] == "parents__name"
    # order is not conserved
    assert df["parents__name"].iloc[0] == [project_label.name]
    # pass a list
    df = ln.ULabel.filter().df(include=["parents__name", "parents__created_by_id"])
    assert df.columns[1] == "parents__created_by_id"
    assert df["parents__name"].iloc[0] == [project_label.name]
    assert set(df["parents__created_by_id"].iloc[0]) == {current_user_id()}

    # for other models
    feature_names = [f"Feature {i}" for i in range(3)]
    features = [ln.Feature(name=name, dtype=int) for name in feature_names]
    ln.save(features)
    feature_set = ln.FeatureSet(features, name="my feature_set")
    feature_set.save()
    feature_set.features.set(features)

    df = ln.FeatureSet.filter(name="my feature_set").df(include="features__name")
    assert df.columns[0] == "features__name"
    # order is not conserved
    assert set(df["features__name"].iloc[0]) == set(feature_names)
    # pass a list
    df = ln.FeatureSet.filter(name="my feature_set").df(
        include=["features__name", "features__created_by_id"]
    )
    assert df.columns[1] == "features__created_by_id"
    assert set(df["features__name"].iloc[0]) == set(feature_names)
    assert set(df["features__created_by_id"].iloc[0]) == {current_user_id()}

    # inner join parents on features
    df = ln.ULabel.filter().df(
        include=["features__name", "features__created_by_id"], join="inner"
    )
    assert set(df["features__name"].iloc[0]) == set(feature_names)
    assert set(df["features__created_by_id"].iloc[0]) == {current_user_id()}
    assert df["parents__name"].iloc[0] == [project_label.name]
    assert set(df["parents__created_by_id"].iloc[0]) == {current_user_id()}

    # outer join parents on features
    df = ln.ULabel.filter().df(
        include=["features__name", "features__created_by_id"], join="outer"
    )
    assert set(df["features__name"].iloc[0]) == set(feature_names)
    assert set(df["features__created_by_id"].iloc[0]) == {current_user_id()}
    assert df["parents__name"].iloc[0] == [project_label.name]
    assert set(df["parents__created_by_id"].iloc[0]) == {current_user_id()}
    print(df)

    # raise error for non many-to-many
    df = ln.ULabel.filter(name="Project 0").df(include="created_by__name")
    assert df["created_by__name"].iloc[0] == "Test User1"

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
    with pytest.raises(NoResultFound):
        qs.one()
    qs = bt.PublicSource.filter().all()
    with pytest.raises(MultipleResultsFound):
        qs.one()
    with pytest.raises(MultipleResultsFound):
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
