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

    dobject = ln.schema.core.dobject(
        id="testid", name="test_file", suffix=".csv", size=1.2
    )
    jupynb = ln.schema.core.jupynb(id="83jf", v="1", name="test")
    ln.db.insert.jupynb(**jupynb.dict())
    dtransform = ln.schema.core.dtransform(jupynb_id=jupynb.id, jupynb_v=jupynb)

    ln.db.insert.dtransform(**dtransform.dict())
    ln.db.insert.dobject_from_dtransform(dobject=dobject, dtransform_id=dtransform.id)

    for entity in ln.schema.list_entities():
        getattr(ln.db.select, entity)().df()

    ln.schema._core.get_db_metadata_as_dict()
    table_object = ln.schema._core.get_table_object("dobject")
    ln.schema._core.get_table_metadata_as_dict(table_object)

    (Path(storage) / "mydata-test-db.lndb").unlink()
    # Note that this merely removes database file but doesn't clean out the instance_settings file!  # noqa
    # Hence, we need to also clean that out:
    from lndb_setup._settings_store import current_instance_settings_file

    current_instance_settings_file.unlink()


if __name__ == "__main__":
    test_create_to_load()
