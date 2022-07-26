from pathlib import Path

from lndb_setup import init_instance, log_in_user

import lamindb as lndb


def test_create_to_load():
    storage = Path.home() / "mydata"
    log_in_user(
        email="raspbear@gmx.de",
        secret="MmR4YuQEyb0yxu7dAwJZTjLzR1Az2lN4Q4IduDlO",
    )
    init_instance(storage=storage)
    lndb.dev.db.insert.dobject(
        name="test_file",
        file_suffix=".csv",
        jupynb_id="83jf",
        jupynb_v="1",
        jupynb_name="test",
        jupynb_type="other",
    )
    for entity in lndb.track.schema.entities:
        lndb.do.load(entity)

    (storage / "mydata.lndb").unlink()


if __name__ == "__main__":
    test_create_to_load()
