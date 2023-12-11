import lamindb as ln


def test_file_visibility():
    # create a file with default visibility
    with open("./test-visibility.txt", "w") as f:
        f.write("visibility")
    artifact = ln.Artifact("./test-visibility.txt", description="test-visibility")
    assert artifact.visibility == 1
    artifact.save()

    # create a dataset from file
    dataset = ln.Dataset(artifact, name="test-visibility")
    dataset.save()

    # delete a dataset will put both dataset and linked file in trash
    dataset.delete()
    assert dataset.artifact.visibility == -1
    result = ln.Dataset.filter(name="test-visibility").all()
    assert len(result) == 0
    result = ln.Dataset.filter(name="test-visibility", visibility=1).all()
    assert len(result) == 0
    result = ln.Dataset.filter(name="test-visibility", visibility=None).all()
    assert len(result) == 1

    # restore
    dataset.restore()
    assert dataset.visibility == 1
    assert dataset.artifact.visibility == 1

    # permanent delete
    dataset.delete(permanent=True)
    result = ln.Artifact.filter(description="test-visibility", visibility=None).all()
    # also permanently deleted linked file
    assert len(result) == 0
