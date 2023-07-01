import lamindb as ln


def test_feature_set_from_values():
    features = ["feat1", "feat2"]
    feature_set = ln.FeatureSet.from_values(features)
    id = feature_set.id
    assert feature_set._state.adding
    assert feature_set.type == "core.Feature"
    assert feature_set.field == "name"
    feature_set.save()
    # test that the feature_set is retrieved from the database
    # in case it already exists
    feature_set = ln.FeatureSet.from_values(features)
    assert not feature_set._state.adding
    assert id == feature_set.id
    feature_set.delete()


def test_feature_set_from_records():
    features = ln.Feature.from_values(["feat1", "feat2"], "name")
    feature_set = ln.FeatureSet(features)
    id = feature_set.id
    assert feature_set._state.adding
    assert feature_set.type == "core.Feature"
    assert feature_set.field == "id"
    feature_set.save()
    # test that the feature_set is retrieved from the database
    # in case it already exists
    feature_set = ln.FeatureSet(features)
    assert not feature_set._state.adding
    assert id == feature_set.id
    feature_set.delete()
