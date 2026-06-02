import lamindb as ln


def test_transfer():
    ln.connect("testdb1")
    ln.Artifact("README.md", key="README.md").save()
    ln.connect("testdb2")
    db1 = ln.DB(f"{ln.setup.settings.user.handle}/testdb1")
    artifact = db1.Artifact.get(key="README.md")
    artifact.save()
    assert artifact.key == "README.md"
    assert artifact.run is not None
    assert artifact.run.finished_at is not None
    assert artifact.run.status == "completed"
    assert artifact.run.started_at is not None
    assert ln.setup.settings.storage.root_as_str.endswith("testdb2")
    assert artifact.storage.root.endswith("testdb1")
