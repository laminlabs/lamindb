import lamindb as ln


def test_query():
    ln.Artifact.connect("laminlabs/lamindata").filter()
