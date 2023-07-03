import pytest

import lamindb as ln


def test_df():
    tag_names = [f"Tag {i}" for i in range(3)]
    projects = [ln.Project(name=f"Project {i}") for i in range(3)]
    ln.save(projects)
    tags = [ln.Tag(name=name) for name in tag_names]
    ln.save(tags)
    for project in projects:
        project.tags.set(tags)
    ln.save(projects)
    df = ln.Project.select().df(include="tags__name")
    assert df.columns[0] == "tags__name"
    # order is not conserved
    assert set(df["tags__name"][0]) == set(tag_names)
    # pass a list
    df = ln.Project.select().df(include=["tags__name", "tags__created_by_id"])
    assert df.columns[1] == "tags__created_by_id"
    assert set(df["tags__name"][0]) == set(tag_names)
    assert set(df["tags__created_by_id"][0]) == set([ln.setup.settings.user.id])

    # raise error for non many-to-many
    with pytest.raises(ValueError):
        ln.Project.select().df(include="external_id")
