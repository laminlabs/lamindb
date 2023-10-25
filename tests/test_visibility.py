import lamindb as ln


def test_file_visibility():
    # create a file with default visibility
    with open("./test-visibility.txt", "w") as f:
        f.write("visibility")
    file = ln.File("./test-visibility.txt", description="test-visibility")
    assert file.visibility == 0
    file.save()

    # filtering
    result = ln.File.filter(description="test-visibility").all()
    assert len(result) == 1

    # file is hidden
    file.visibility = 1
    file.save()
    result = ln.File.filter(description="test-visibility").all()
    assert len(result) == 0
    result = ln.File.filter(description="test-visibility", visibility=None).all()
    assert len(result) == 1
    result = ln.File.filter(description="test-visibility", visibility="default").all()
    assert len(result) == 0
