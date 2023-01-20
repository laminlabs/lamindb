from pathlib import Path


def test_create_to_load():
    import lndb_setup

    storage_root = "mydata-test-db"
    lndb_setup.login(
        "raspbear@gmx.de",
        password="MmR4YuQEyb0yxu7dAwJZTjLzR1Az2lN4Q4IduDlO",
    )
    lndb_setup.init(storage=storage_root)

    # if importing this at the top of the file lndb_setup will try to
    # create unnecessary tables
    import lamindb as ln

    notebook = ln.schema.Notebook(id="83jf", v="1", name="test")
    ln.add(notebook)
    run = ln.schema.Run(notebook_id=notebook.id, notebook_v=notebook.v)
    ln.add(run)
    ln.select(ln.schema.Storage, root=str(ln.settings.instance.storage.root)).one()

    ln.schema._core.get_db_metadata_as_dict()
    table_object = ln.schema._core.get_table_object("core.dobject")
    ln.schema._core.get_table_metadata_as_dict(table_object)

    (Path(storage_root) / "mydata-test-db.lndb").unlink()
    # Note that this merely removes database file but doesn't clean out the instance_settings file!  # noqa
    # Hence, we need to also clean that out:
    from lndb_setup._settings_store import current_instance_settings_file

    current_instance_settings_file().unlink()


if __name__ == "__main__":
    test_create_to_load()
