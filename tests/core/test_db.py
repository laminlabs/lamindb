import lamindb as ln


def test_create_to_load():
    transform = ln.Transform(version="0", key="test", type="pipeline")
    transform.save()
    run = ln.Run(transform=transform)
    run.save()
    ln.Storage.get(root=str(ln.setup.settings.storage.root))
