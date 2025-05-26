import lamindb as ln


def testbranch_id():
    # create a file with default branch_id
    with open("./testbranch_id.txt", "w") as f:
        f.write("branch_id")
    artifact = ln.Artifact("./testbranch_id.txt", description="testbranch_id").save()
    assert artifact.branch_id == 1

    # create a collection from file
    collection = ln.Collection(artifact, key="testbranch_id").save()

    # delete a collection will put both collection but not linked artifact in trash
    collection.delete()
    assert collection.ordered_artifacts[0].branch_id == 1
    result = ln.Collection.filter(key="testbranch_id").all()
    assert len(result) == 0
    result = ln.Collection.filter(key="testbranch_id", branch_id=1).all()
    assert len(result) == 0
    result = ln.Collection.filter(key="testbranch_id", visibility=1).all()
    assert len(result) == 0
    result = ln.Collection.filter(key="testbranch_id", branch_id=None).all()
    assert len(result) == 1
    result = ln.Collection.filter(key="testbranch_id", visibility=None).all()
    assert len(result) == 1

    # restore
    collection.restore()
    assert collection.branch_id == 1
    assert collection.ordered_artifacts[0].branch_id == 1

    # permanent delete
    collection.delete(permanent=True)
    result = ln.Artifact.filter(description="testbranch_id", branch_id=None).all()
    # also permanently deleted linked file
    assert len(result) == 1
