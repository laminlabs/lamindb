import lamindb as ln


def test_manager_list():
    tag = ln.Tag(name="manager tag")
    tag.save()
    tag_names = [f"Tag {i}" for i in range(3)]
    tags = [ln.Tag(name=name) for name in tag_names]
    ln.save(tags)
    tag.parents.set(tags)
    assert len(tag.parents.list()) == 3
    assert "Tag 1" in tag.parents.list("name")
    tag.delete()
    for tag in tags:
        tag.delete()
