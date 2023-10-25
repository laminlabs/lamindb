import lamindb as ln


def test_file_visibility():
    # create a file with default visibility
    with open("./test-visibility.txt", "w") as f:
        f.write("visibility")
    file = ln.File("./test-visibility.txt", description="test-visibility")
    assert file.visibility == 0
    file.save()

    # create a dataset from file
    dataset = ln.Dataset(file, name="test-visibility")
    dataset.save()

    # delete a dataset will put both dataset and linked file in trash
    dataset.delete()
    assert dataset.file.visibility == 2
    result = ln.Dataset.filter(description="test-visibility").all()
    assert len(result) == 0
    result = ln.Dataset.filter(
        description="test-visibility", visibility="default"
    ).all()
    assert len(result) == 0
    result = ln.Dataset.filter(description="test-visibility", visibility=None).all()
    assert len(result) == 1

    # delete a file
    result = ln.File.filter(description="test-visibility").all()
    assert len(result) == 0
    result = ln.File.filter(description="test-visibility", visibility=None).all()
    assert len(result) == 1

    # delete a dataset from trash
    dataset.delete(force=True)
    result = ln.File.filter(description="test-visibility", visibility=None).all()
    assert len(result) == 0
