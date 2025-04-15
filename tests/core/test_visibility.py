import lamindb as ln


def test_branch_code():
    # create a file with default _branch_code
    with open("./test_branch_code.txt", "w") as f:
        f.write("_branch_code")
    artifact = ln.Artifact(
        "./test_branch_code.txt", description="test_branch_code"
    ).save()
    assert artifact._branch_code == 1

    # create a collection from file
    collection = ln.Collection(artifact, key="test_branch_code").save()

    # delete a collection will put both collection but not linked artifact in trash
    collection.delete()
    assert collection.ordered_artifacts[0]._branch_code == 1
    result = ln.Collection.filter(key="test_branch_code").all()
    assert len(result) == 0
    result = ln.Collection.filter(key="test_branch_code", _branch_code=1).all()
    assert len(result) == 0
    result = ln.Collection.filter(key="test_branch_code", visibility=1).all()
    assert len(result) == 0
    result = ln.Collection.filter(key="test_branch_code", _branch_code=None).all()
    assert len(result) == 1
    result = ln.Collection.filter(key="test_branch_code", visibility=None).all()
    assert len(result) == 1

    # restore
    collection.restore()
    assert collection._branch_code == 1
    assert collection.ordered_artifacts[0]._branch_code == 1

    # permanent delete
    collection.delete(permanent=True)
    result = ln.Artifact.filter(description="test_branch_code", _branch_code=None).all()
    # also permanently deleted linked file
    assert len(result) == 1
