import lamindb as ln


def test_create_to_load():
    transform = ln.Transform(version="0", name="test", type="pipeline")
    ln.add(transform)
    run = ln.Run(transform=transform)
    ln.add(run)
    ln.select(ln.Storage, root=str(ln.setup.settings.storage.root)).one()

    ln.schema._core.get_db_metadata_as_dict()
    table_object = ln.schema._core.get_table_object("lnschema_core_file")
    ln.schema._core.get_table_metadata_as_dict(table_object)


if __name__ == "__main__":
    test_create_to_load()
