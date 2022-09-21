from pathlib import Path


def test_create_to_load():
    import lndb_setup

    storage = "mydata-test-db"
    lndb_setup.login(
        "raspbear@gmx.de",
        password="MmR4YuQEyb0yxu7dAwJZTjLzR1Az2lN4Q4IduDlO",
    )
    lndb_setup.init(storage=storage)

    # if importing this at the top of the file lndb_setup will try to
    # create unnecessary tables
    import lamindb as ln

    ln.db.insert.dobject_from_jupynb(
        name="test_file",
        file_suffix=".csv",
        jupynb_id="83jf",
        jupynb_v="1",
        jupynb_name="test",
    )
    for entity in ln.schema.list_entities():
        getattr(ln.db.query, entity)(as_df=True).all()

    (Path(storage) / "mydata-test-db.lndb").unlink()
    # Note that this merely removes database file but doesn't clean out the instance_settings file!  # noqa
    # Hence, we need to also clean that out:
    from lndb_setup._settings_store import current_instance_settings_file

    current_instance_settings_file.unlink()


if __name__ == "__main__":
    test_create_to_load()
