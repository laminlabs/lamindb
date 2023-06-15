import lamindb as ln


def test_create_to_load():
    transform = ln.Transform(version="0", name="test", type="pipeline")
    ln.save(transform)
    run = ln.Run(transform=transform)
    ln.save(run)
    ln.select(ln.Storage, root=str(ln.setup.settings.storage.root)).one()
