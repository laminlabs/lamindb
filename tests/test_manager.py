import lamindb as ln


def test_manager_list():
    project = ln.Project(name="manager project")
    project.save()
    tag_names = [f"Tag {i}" for i in range(3)]
    tags = [ln.Tag(name=name) for name in tag_names]
    ln.save(tags)
    project.tags.set(tags)
    assert len(project.tags.list()) == 3
    assert "Tag 1" in project.tags.list("name")
