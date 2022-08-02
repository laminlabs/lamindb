from pathlib import Path

import lndb_setup

import lamindb as db


def test_create_to_load():
    storage = "mydata-test-db"
    lndb_setup.login(
        "raspbear@gmx.de",
        password="MmR4YuQEyb0yxu7dAwJZTjLzR1Az2lN4Q4IduDlO",
    )
    lndb_setup.init(storage=storage)
    db.dev.db.insert.dobject_from_jupynb(
        name="test_file",
        file_suffix=".csv",
        jupynb_id="83jf",
        jupynb_v="1",
        jupynb_name="test",
    )
    for entity in db.track.schema.entities:
        db.do.load.entity(entity)

    (Path(storage) / "mydata-test-db.lndb").unlink()
    # Note that this merely removes database file but doesn't clean out the instance_settings file!  # noqa
    # Hence, we need to also clean that out:
    from lndb_setup._settings_store import current_instance_settings_file

    current_instance_settings_file.unlink()


if __name__ == "__main__":
    test_create_to_load()
