import lamindb as ln


def test_file_visibility():
    # create a file with default visibility
    with open("./test-visibility.txt", "w") as f:
        f.write("visibility")
    artifact = ln.Artifact("./test-visibility.txt", description="test-visibility")
    assert artifact.visibility == 1
    artifact.save()

    # create a collection from file
    collection = ln.Collection(artifact, name="test-visibility")
    collection.save()

    # delete a collection will put both collection but not linked artifact in trash
    collection.delete()
    assert collection.ordered_artifacts[0].visibility == 1
    result = ln.Collection.filter(name="test-visibility").all()
    assert len(result) == 0
    result = ln.Collection.filter(name="test-visibility", visibility=1).all()
    assert len(result) == 0
    result = ln.Collection.filter(name="test-visibility", visibility=None).all()
    assert len(result) == 1

    # restore
    collection.restore()
    assert collection.visibility == 1
    assert collection.ordered_artifacts[0].visibility == 1

    # permanent delete
    collection.delete(permanent=True)
    result = ln.Artifact.filter(description="test-visibility", visibility=None).all()
    # also permanently deleted linked file
    assert len(result) == 1
