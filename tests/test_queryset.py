import pytest

import lamindb as ln


def test_df():
    project_tag = ln.Tag(name="Project")
    project_tag.save()
    project_names = [f"Project {i}" for i in range(3)]
    tags = [ln.Tag(name=name) for name in project_names]
    ln.save(tags)
    for tag in tags:
        tag.parents.add(project_tag)
    df = ln.Tag.select().df(include="parents__name")
    assert df.columns[0] == "parents__name"
    # order is not conserved
    assert df["parents__name"][0] == [project_tag.name]
    # pass a list
    df = ln.Tag.select().df(include=["parents__name", "parents__created_by_id"])
    assert df.columns[1] == "parents__created_by_id"
    assert df["parents__name"][0] == [project_tag.name]
    assert set(df["parents__created_by_id"][0]) == set([ln.setup.settings.user.id])

    # raise error for non many-to-many
    with pytest.raises(ValueError):
        ln.Tag.select().df(include="name")

    # call it from a non-select-derived queryset
    qs = ln.User.objects.all()
    assert qs.df().iloc[0]["handle"] == "testuser1"


def test_one_first():
    qs = ln.User.objects.all()
    assert qs.one().handle == "testuser1"
    assert qs.first().handle == "testuser1"
    assert qs.one_or_none().handle == "testuser1"


def test_search():
    tag_names = [f"Tag {i}" for i in range(3)]
    tags = [ln.Tag(name=name) for name in tag_names]
    ln.save(tags)
    qs = ln.Tag.select(name="Tag 2").all()
    assert qs.search("Tag 1").iloc[0].name == "Tag 2"


def test_lookup():
    qs = ln.User.select(handle="testuser1").all()
    lookup = qs.lookup(field="handle")
    assert lookup.testuser1.handle == "testuser1"


def test_inspect():
    qs = ln.User.select(handle="testuser1").all()
    assert qs.inspect(["user1", "user2"], "name")["mapped"] == []


def test_map_synonyms():
    qs = ln.User.select(handle="testuser1").all()
    assert qs.map_synonyms(["user1", "user2"]) == ["user1", "user2"]
