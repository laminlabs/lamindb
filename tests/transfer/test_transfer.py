import lamindb as ln


def test_transfer():
    ln.connect("testdb1")
    ln.Artifact("README.md", key="README.md").save()
    ln.connect("testdb2")
    db1 = ln.DB(f"{ln.setup.settings.user.handle}/testdb1")
    artifact = db1.Artifact.get(key="README.md")
    artifact.save()
    assert artifact.key == "README.md"
    assert ln.setup.settings.storage.root_as_str.endswith("testdb2")
    assert artifact.storage.root.endswith("testdb1")
