import lamindb as ln


def test_file__branch_code():
    # create a file with default _branch_code
    with open("./test-_branch_code.txt", "w") as f:
        f.write("_branch_code")
    artifact = ln.Artifact("./test-_branch_code.txt", description="test-_branch_code")
    assert artifact._branch_code == 1
    artifact.save()

    # create a collection from file
    collection = ln.Collection(artifact, name="test-_branch_code")
    collection.save()

    # delete a collection will put both collection but not linked artifact in trash
    collection.delete()
    assert collection.ordered_artifacts[0]._branch_code == 1
    result = ln.Collection.filter(name="test-_branch_code").all()
    assert len(result) == 0
    result = ln.Collection.filter(name="test-_branch_code", _branch_code=1).all()
    assert len(result) == 0
    result = ln.Collection.filter(name="test-_branch_code", _branch_code=None).all()
    assert len(result) == 1

    # restore
    collection.restore()
    assert collection._branch_code == 1
    assert collection.ordered_artifacts[0]._branch_code == 1

    # permanent delete
    collection.delete(permanent=True)
    result = ln.Artifact.filter(
        description="test-_branch_code", _branch_code=None
    ).all()
    # also permanently deleted linked file
    assert len(result) == 1
